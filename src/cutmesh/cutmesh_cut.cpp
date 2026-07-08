// File for mesh cutting
#include "cutmesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <utility>

// File for mesh cutting operations

namespace Mesh {
    // Recursively computes straightest geodesic from a starting point given a starting direction, inserting new vertices as needed
    int cutmesh::embedCurves() {
        if (dCN == nullptr || M == nullptr) {
            return -1;
        }
        std::map<int, int> dCNVtoV;
        const std::vector<DCurvenet::Vert>& dCN_V = dCN->V;
        std::vector<Vert> proj_V(dCN_V.size());
        // 1. First, create a list of all projections (parallel)
        #pragma omp parallel for
        for (int v = 0; v < dCN_V.size(); v++) {
            Eigen::Vector3d pos = dCN_V[v].pos;
            Eigen::Vector3d proj;
            Eigen::Vector3d n;
            int elIdx;
            int elType = M->computeVProjection(pos, proj, elIdx);
            if (elType == 0) {  // Vertex
                n = M->V[elIdx].n;
            } else if (elType == 1) {   // Edge
                n = M->E[elIdx].n;
            } else {    // Face
                n = M->F[elIdx].n;
            }
            proj_V[v] = createVertex(proj, n, 1, -1, elType, elIdx, pos-proj);
        }
        
        // Now, add all vertices into the cutmesh (SEQUENTIAL)
        std::map<int, std::vector<std::pair<double, int>>> edgeMap;
        for (int v = 0; v < proj_V.size(); v++) {
            // First, identify what kind of vertex should be inserted
            int elType = proj_V[v].mesh_elType;
            int elIdx = proj_V[v].mesh_elIdx;
            // If it landed on a vertex, then modify the existing vertex
            if (elType == 0) {
                V[elIdx] = proj_V[v];
                dCNVtoV[v] = elIdx;
            } else if (elType == 1) {
                // If it landed on an edge, we need to figure out where exacty to split the edge
                int new_v = insertVertex(proj_V[v]);
                // edgeMap[elIdx].push_back(new_v);
                dCNVtoV[v] = new_v;
                // Check where to split
                int insert_index = 0;
                int e_insert = elIdx;
                Eigen::Vector3d v0 = M->V[HE[E[elIdx].he].dest].pos;
                Eigen::Vector3d v1 = M->V[HE[HE[E[elIdx].he].twin].dest].pos;
                double t = (V[new_v].pos - v0).norm()/(v1 - v0).norm();
                if (edgeMap.find(elIdx) != edgeMap.end()) { // Exists in map, we need to compute t vals etc.
                    std::vector<std::pair<double, int>> localSplits = edgeMap[elIdx];
                    // Find where to insert the new vertex
                    for (int j = 0; j < localSplits.size(); j++) {
                        if (localSplits[j].first <= t) {
                            insert_index++;
                        }
                    }
                    int e_insert = 0;
                    // Insert into mesh via edge split
                    if (insert_index != 0) {
                        e_insert = localSplits[insert_index - 1].second;
                    }
                }
                int new_e = splitEdge(e_insert, new_v);
                edgeMap[elIdx].insert(edgeMap[elIdx].begin() + insert_index, std::make_pair(t, new_e));
            } else {
                // If it landed on a face, then simply insert the vertex
                int new_v = insertVertex(proj_V[v]);
                dCNVtoV[v] = new_v;
            }
        }

        // Trace the curves (SEQUENTIAL, due to geodesics)
        const std::vector<DCurvenet::Edge>& dCN_E = dCN->E;
        const std::vector<DCurvenet::HalfEdge>& dCN_HE = dCN->HE;
        for (int e = 0; e < dCN_E.size(); e++) {
            // For this edge, connect the two associated vertices in the cutmesh
            int dCN_he0 = dCN_E[e].he;
            int dCN_he1 = dCN_HE[dCN_he0].twin;
            int dCN_v0 = dCN_HE[dCN_he1].dest;
            int dCN_v1 = dCN_HE[dCN_he0].dest;
            int v0 = dCNVtoV[dCN_HE[dCN_he1].dest];
            int v1 = dCNVtoV[dCN_HE[dCN_he0].dest];

            // If edge already exists, modify the existing edge
            if (vertPairToHE.find(std::make_pair(v0, v1)) != vertPairToHE.end()) {
                int he0 = vertPairToHE[std::make_pair(v0, v1)];
                int he1 = HE[he0].twin;
                HE[he0].dCN_idx = dCN_he0;
                HE[he1].dCN_idx = dCN_he1;
                continue;
            }
            // Classify based on properties
            std::vector<Vert> traceVerts;
            traceVerts.push_back(V[v0]);
            Eigen::Vector3d direc = (V[v1].pos - V[v0].pos).normalized();
            int success = M->traceGeodesic(V[v0], V[v1], direc, 
                                           V[v0].mesh_elType, V[v0].mesh_elIdx, 
                                           traceVerts);
            // This is a likely spot for failure, so flag it
            if (success == -1) {
                return -1;
            }
            traceVerts.push_back(V[v1]);

            // Insert each vert into the cutmesh
            std::vector<int> traceList;
            traceList.push_back(v0);
            for (int v_idx = 1; v_idx < traceVerts.size() - 1; v_idx++) {
                traceVerts[v_idx].label = 2;
                // TODO: Is this a safe thing to do?
                // i.e., the projVector taken using the midpoint of the associated halfedge, thus being deterministic without 
                // Requiring us to subdivide the dCN edge OR change the projection estimate formula (sort of)
                traceVerts[v_idx].projVector = ((dCN_V[dCN_v1].pos + dCN_V[dCN_v0].pos) / 2) - traceVerts[v_idx].pos;
                // First, compute an estimated curvenet position so we can take the difference
                if (traceVerts[v_idx].mesh_elType == 0) {  // Check if we are on a vertex
                    V[traceVerts[v_idx].mesh_elIdx] = traceVerts[v_idx];
                    traceList.push_back(traceVerts[v_idx].mesh_elIdx);
                } else { // We must be on an edge
                    int new_v = insertVertex(traceVerts[v_idx]);
                    // Check if the edge was already split. If so, find where to split it.
                    int insert_index = 0;
                    int orig_e = V[new_v].mesh_elIdx;   // Edge from the original mesh
                    int e_insert = orig_e;      // Edge in cut-mesh to insert at
                    Eigen::Vector3d v0 = M->V[HE[E[e_insert].he].dest].pos;
                    Eigen::Vector3d v1 = M->V[HE[HE[E[e_insert].he].twin].dest].pos;
                    double t = (V[new_v].pos - v0).norm()/(v1 - v0).norm();
                    if (edgeMap.find(orig_e) != edgeMap.end()) { // Exists in map, we need to compute t vals etc.
                        std::vector<std::pair<double, int>> localSplits = edgeMap[orig_e];
                        // Find where to insert the new vertex
                        for (int j = 0; j < localSplits.size(); j++) {
                            if (localSplits[j].first <= t) {
                                insert_index++;
                            }
                        }
                        int e_insert = 0;
                        // Insert into mesh via edge split
                        if (insert_index != 0) {
                            e_insert = localSplits[insert_index - 1].second;
                        }
                    }
                    int new_e = splitEdge(e_insert, new_v);
                    edgeMap[orig_e].insert(edgeMap[orig_e].begin() + insert_index, std::make_pair(t, new_e));
                    traceList.push_back(new_v);
                }
            }
            traceList.push_back(v1);

            // Connect the inserted vertices
            for (int v_idx = 1; v_idx < traceVerts.size(); v_idx++) {
                int e = insertEdge(traceList[v_idx - 1], traceList[v_idx], dCN_he0, dCN_he1);
                if (e == -1) {  // Edge already exists
                    int he0 = vertPairToHE[{v0, v1}];
                    int he1 = vertPairToHE[{v1, v0}];
                    HE[he0].dCN_idx = dCN_he0;
                    HE[he1].dCN_idx = dCN_he1;
                }
            }
        }

        // Compute faces
        sortHalfEdges();
        resetFaces();

        // "Remove" obsolete curves by deactivating them
        deactivateIsolatedCuts();

        return 1;
    }

    // Properly orient halfedges based on local tangent planes
    int cutmesh::sortHalfEdges() {
        // Determine incoming halfedges for all vertices
        std::vector<std::vector<int>> outHE(V.size());
        for (int he = 0; he < HE.size(); he++) {
            if (!HE[he].active) {
                continue;
            }
            outHE[HE[HE[he].dest].twin].push_back(he);
        }

        // Iterate over vertices and sort their halfedges based on projection to the tangent plane
        for (int v = 0; v < V.size(); v++) {
            if (!V[v].active) {
                continue;
            }
            std::vector<int> sortedHE = outHE[v];
            if (sortedHE.size() == 1) { // Handle endpoints
                HE[sortedHE[0]].prev = HE[sortedHE[0]].twin;
                HE[HE[sortedHE[0]].twin].next = sortedHE[0];
                continue;
            }
            // Project all onto local tangent plane and sort
            std::vector<double> projAngles(sortedHE.size());
            Eigen::Vector3d t1, t2;
            Utils::buildPlaneBasis(V[v].n, t1, t2);
            for (int he = 0; he < sortedHE.size(); he++) {
                V[HE[sortedHE[he]].dest].pos;
                double theta;
                Utils::directionAngleInPlane(V[v].pos, V[HE[sortedHE[he]].dest].pos, V[v].n, t1, t2, theta);
                projAngles[he] = theta;
            }
            Utils::doubleListIdxSort(projAngles, sortedHE);
            
            // Rewire based on sorted order
            for (int i = 0; i < sortedHE.size(); i++) {
                int he = sortedHE[i];
                int he_p1 = sortedHE[(i+1)%sortedHE.size()];
                
                HE[he].prev = HE[he_p1].twin;
                HE[he_p1].next = he;
            }
        }
        return 1;
    }

    // Reset every face based on the halfedges
    int cutmesh::resetFaces() {
        // Throw away all old faces
        F.clear();
        // Iterate over halfedges and find every loop
        std::vector<bool> seenHE(HE.size());
        for (int he = 0; he < HE.size(); he++) {
            // Ignore inactive and boundary
            if (seenHE[he] || !HE[he].active || HE[he].boundary) {
                seenHE[he] = true;
                continue;
            }
            int curr_he = he;
            int counter = 1;
            std::vector<int> heLoop;
            // Loop until we get the face
            // Push back the face
            int f = F.size();
            do {
                heLoop.push_back(curr_he);
                seenHE[he] = true;
                HE[curr_he].face = f;
                curr_he = HE[curr_he].next;
                counter++;
            } while (curr_he != he || counter >= HE.size());
            F.emplace_back();
            F[f].he = he;
        }
        // Compute areas and normals for all created faces
        computeFNormalsAreas();
        return 1;
    }

    // "Clean up" curves by deactivating any edges, halfedges, faces, and verts that are not useful/fully enclosed (see paper)
    int cutmesh::deactivateIsolatedCuts() {
        // Iterate over faces and deactivate if all vertices are projected onto mesh
        for (int f = 0; f < F.size(); f++) {
            std::vector<int> adjHE = faceAdjHalfEdges(f);
            bool deactivate = true;
            for (int he = 0; he < adjHE.size(); he++) {
                int v = HE[adjHE[he]].dest;
                if ((V[v].mesh_elType != 2)) {
                    // At least one vertex in the face is either from the original mesh, landed on a mesh vertex, or split a mesh edge
                    deactivate = false;
                }
            }

            // deactivate all vertices, edges and faces in the loop
            if (deactivate) {
                for (int he = 0; he < adjHE.size(); he++) {
                    HE[adjHE[he]].active = false;
                    HE[HE[adjHE[he]].twin].active = false;
                    V[HE[adjHE[he]].dest].active = false;
                    E[HE[adjHE[he]].edge].active = false;
                }
                F[f].active = false;
            }
        }

        // Loose curves that form no loops.
        // Iterate over halfedges and identify which are without a face and not boundary
        for (int he = 0; he < HE.size(); he++) {
            // NOTE: face index of -1 indicates boundary, meaning that an additional boundary flag set to false is a contradiction
            if (HE[he].face == -1 && !HE[he].boundary) {
                // Deactivate the associated vertices and edge
                V[HE[he].dest].active = false;
                E[HE[he].edge].active = false;
                HE[he].active = false;
            }
        }

        // Loose verts (for sanity)
        for (int v = 0; v < V.size(); v++) {
            if (V[v].active && (vertAdjHEs(v).size() == 0)) {
                V[v].active = false;
            }
        }

        return 1;
    }

    // Cuts mesh by "unzipping" along curve vertices
    int cutmesh::cutMesh() {
        // Iterate over each dCNVert
        for (int v = 0; v < V.size(); v++) {
            if (V[v].label != 0) {  // Only process dCN vertices
                // Gather all adjacent halfedges
                std::vector<int> adjHE = vertAdjHEs(v);
                // Identify each dCN he in the local list
                std::vector<int> localHE_split;
                int num_dCN = 0;
                for (int he = 0; he < adjHE.size(); he++) {
                    // Track local halfedges that come from the dCN or are at boundaries
                    if (HE[adjHE[he]].dCN_idx != -1) {
                        localHE_split.push_back(he);
                        num_dCN++;
                    } else if (HE[adjHE[he]].boundary) {
                        localHE_split.push_back(he);
                    }
                }
                // If there is only one outgoing, then skip
                if (num_dCN == 0) {
                    continue;
                } else if (num_dCN == 1 && !vertIsBoundary(v)) {
                    // If vert is not on the boundary and only has one adjacent cut edge, then don't split
                    V[v].corner_idx = HE[localHE_split[0]].twin;
                    continue;
                }

                // Find sequences of dCN corners
                std::vector<std::vector<int>> dCN_corners(localHE_split.size());
                for (int i = 0; i < localHE_split.size(); i++) {
                    int globalStart = adjHE[localHE_split[i]];
                    int globalNext = adjHE[localHE_split[(i+1)%localHE_split.size()]];
                    std::vector<int> he_corner;
                    int globalCurrent = globalStart;
                    do {
                        he_corner.push_back(globalCurrent);
                        he_corner.push_back(HE[globalCurrent].prev);
                        globalCurrent = HE[HE[globalCurrent].prev].twin;
                    } while (globalCurrent != globalNext);
                    dCN_corners[i] = he_corner;
                }
                // Check for existing boundary halfedges
                for (int c = 0; c < dCN_corners.size(); c++) {
                    if (HE[dCN_corners[c][0]].boundary) {   // If the first one is a boundary
                        
                    }
                }


                // For each dCN corner, split the vertex and rewire the halfedge destination vertices accordingly
                bool first = true;
                for (int i = 0; i < localHE_split.size(); i++) {
                    int localStart = localHE_split[i];
                    int localNext = localHE_split[(i+1)%localHE_split.size()];
                    // Incoming next
                    int he1 = HE[adjHE[localNext]].twin;
                    // Outgoing start
                    int he0 = adjHE[localStart];
                    // Skip duplication if both are on the boundary
                    if (HE[he0].boundary && HE[he1].boundary) {
                        continue;
                    }
                    int v_current;
                    if (!first) {
                        v_current = insertVertex(V[v].pos, V[v].n, V[v].label, -1, V[v].mesh_elType, V[v].mesh_elIdx, V[v].projVector);
                    } else {
                        v_current = v;
                        first = false;
                    }
                    // Rewire corner index
                    V[v_current].corner_idx = he1;
                    V[v_current].he = he0;
                    int he_current = localStart;
                    // Rewire the halfedges to point to the new vert
                    do {
                        int adj_v = HE[HE[he_current].twin].dest;
                        HE[HE[he_current].prev].dest = v_current;
                        // Reset the vertPairToHE
                        vertPairToHE.erase({v, });
                        vertPairToHE.erase({});
                        vertPairToHE[{u, v_current}] = he0_idx;
                        vertPairToHE[{v_current, }] = he0_idx;
                        he_current = HE[HE[he_current].prev].twin;
                    } while(he_current != localNext);
                }
            }
        }

        // Loop over edges and "cut" them apart if necessary
        // TODO: Potentially fill these in with new boundary halfedges rather than nothing
        for (int e = 0; e < E.size(); e++) {
            int he0 = E[e].he;
            int he1 = HE[he0].twin;
            if (HE[he0].dCN_idx == -1) { // Process only those from dCN
                continue;
            }
            // Separate the halfedges
            HE[he0].twin = -1;
            HE[he1].twin = -1;
        }

        return 1;
    }

}   // namespace Mesh