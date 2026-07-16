// File for mesh cutting
#include "cutmesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <utility>
#include <iostream>

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
            vertProjData projData = M->computeVProjection(pos, proj);
            if (projData.elType == 0) {  // Vertex
                n = M->V[projData.elIdx].n;
            } else if (projData.elType == 1) {   // Edge
                n = M->E[projData.elIdx].n;
            } else {    // Face
                n = M->F[projData.elIdx].n;
            }
            proj_V[v] = createVertex(proj, n, 1, -1, projData.elType, projData.elIdx, pos-proj);
        }
        
        // Now, add all vertices into the cutmesh (SEQUENTIAL)
        std::map<int, std::vector<std::pair<double, int>>> edgeMap;
        for (int v = 0; v < proj_V.size(); v++) {
            // First, identify what kind of vertex should be inserted
            int elType = proj_V[v].projData.elType;
            int elIdx = proj_V[v].projData.elIdx;
            // If it landed on a vertex, then modify the existing vertex
            if (elType == 0) {
                // If we landed on an existing cut-vertex, then things are bad!
                if (V[elIdx].label != 0) {
                    return -1;
                }
                V[elIdx].label = 1;
                V[elIdx].projData = proj_V[v].projData;
                V[elIdx].defData.projVector = proj_V[v].defData.projVector;
                dCNVtoV[v] = elIdx;
            } else if (elType == 1) {
                // If it landed on an edge, we need to figure out where exacty to split the edge
                int new_v = insertVertex(proj_V[v]);
                dCNVtoV[v] = new_v;
                // Check where to split
                int insert_index = 0;
                int e_insert = elIdx;
                int he = M->E[elIdx].he;
                int a = M->HE[M->HE[he].twin].dest;
                int b = M->HE[he].dest;

                Eigen::Vector3d p0 = M->V[a].pos;
                Eigen::Vector3d p1 = M->V[b].pos;
                double t = (V[new_v].pos - p0).norm() / (p1 - p0).norm();
                if (edgeMap.find(elIdx) != edgeMap.end()) { // Exists in map, we need to compute t vals etc.
                    const std::vector<std::pair<double, int>>& localSplits = edgeMap[elIdx];
                    // Find where to insert the new vertex
                    for (int j = 0; j < localSplits.size(); j++) {
                        if (std::abs(localSplits[j].first - t) <= 1e-12) {    // Landed on the same point as a different vertex
                            return -1;
                        } else if (localSplits[j].first < t) {
                            insert_index++;
                        }
                    }
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
            int depth = 0;
            Eigen::Vector3d direc = (V[v1].pos - V[v0].pos).normalized();
            //Eigen::Vector3d direc = (dCN_V[dCN_v1].pos - dCN_V[dCN_v0].pos).normalized();
            std::cout << "Starting trace from vert " << v0 << " (proj type: " << V[v0].projData.elType << ") to " << v1 << " (proj type: " << V[v1].projData.elType << ")" << std::endl;
            std::cout << "Initial direction: " << direc[0] << ", "
                                        << direc[1] << ", "
                                        << direc[2] << std::endl;
            int success = M->traceGeodesic(V[v0], V[v1], direc, 
                                           V[v0].projData,
                                           traceVerts, depth);
            // This is a likely spot for failure, so flag it
            if (success != 1) {
                std::cout << "Trace failed. Trying opposite direction." << std::endl;
                // Check the opposite direction to catch initial directional error
                // TODO: How necessary is this?
                depth = 0;
                traceVerts.clear();
                traceVerts.push_back(V[v0]);
                int success_opposite = M->traceGeodesic(V[v0], V[v1], -1 * direc, 
                                           V[v0].projData,
                                           traceVerts, depth);
                if (success_opposite != 1) {
                    std::cout << "Trace failed again. Quitting." << std::endl;
                    return -1;
                }
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
                traceVerts[v_idx].defData.projVector = ((dCN_V[dCN_v1].pos + dCN_V[dCN_v0].pos) / 2) - traceVerts[v_idx].pos;
                // First, compute an estimated curvenet position so we can take the difference
                if (traceVerts[v_idx].projData.elType == 0) {  // Check if we are on a vertex
                    if (V[traceVerts[v_idx].projData.elIdx].label != 0) { // If we hit a vertex that is already assigned, then quit
                        return -1;
                    }
                    V[traceVerts[v_idx].projData.elIdx].label = 2;
                    V[traceVerts[v_idx].projData.elIdx].projData = traceVerts[v_idx].projData;
                    V[traceVerts[v_idx].projData.elIdx].defData = traceVerts[v_idx].defData;
                    traceList.push_back(traceVerts[v_idx].projData.elIdx);
                } else if (traceVerts[v_idx].projData.elType == 1) { // We must be on an edge
                    int new_v = insertVertex(traceVerts[v_idx]);
                    // Check if the edge was already split. If so, find where to split it.
                    int orig_e = V[new_v].projData.elIdx;   // source mesh edge index
                    int e_insert = orig_e;                  // cutmesh edge to split
                    // Compute t ONLY from the original/source edge
                    int orig_he = M->E[orig_e].he;
                    int a = M->HE[M->HE[orig_he].twin].dest;
                    int b = M->HE[orig_he].dest;

                    Eigen::Vector3d p0 = M->V[a].pos;
                    Eigen::Vector3d p1 = M->V[b].pos;
                    double t = (V[new_v].pos - p0).norm() / (p1 - p0).norm();

                    int insert_index = 0;
                    if (edgeMap.find(orig_e) != edgeMap.end()) { // Exists in map, we need to compute t vals etc.
                        const std::vector<std::pair<double, int>>& localSplits = edgeMap[orig_e];
                        // Find where to insert the new vertex
                        for (int j = 0; j < localSplits.size(); j++) {
                            if (std::abs(localSplits[j].first - t) <= 1e-12) {
                                return -1;
                            } else if (localSplits[j].first < t) {
                                insert_index++;
                            }
                        }
                        // Insert into mesh via edge split
                        if (insert_index != 0) {
                            e_insert = localSplits[insert_index - 1].second;
                        }
                    }
                    int new_e = splitEdge(e_insert, new_v);
                    edgeMap[orig_e].insert(edgeMap[orig_e].begin() + insert_index, std::make_pair(t, new_e));
                    traceList.push_back(new_v);
                } else {    // Something weird happened
                    return -1;
                }
            }
            traceList.push_back(v1);
            // Connect the inserted vertices
            for (int v_idx = 1; v_idx < traceVerts.size(); v_idx++) {
                int v_start = traceList[v_idx - 1];
                int v_next = traceList[v_idx];
                int e = insertEdge(v_start, v_next, dCN_he0, dCN_he1);
                if (e == -1) {  // Edge already exists
                    int he0 = vertPairToHE[{v_start, v_next}];
                    int he1 = vertPairToHE[{v_next, v_start}];
                    HE[he0].dCN_idx = dCN_he0;
                    HE[he1].dCN_idx = dCN_he1;
                }
            }
        }
        std::cout << "Inserted all new vertices and edges" << std::endl;
        // Compute faces
        if (sortHalfEdges() != 1) {
            return -1;
        }
        std::cout << "Halfedges sorted." << std::endl;
        if (resetFaces() != 1) {
            return -1;
        }
        std::cout << "Faces reset." << std::endl;
        // "Remove" obsolete curves by deactivating them
        //deactivateIsolatedCuts();
        //std::cout << "Isolated cuts deactivated." << std::endl;
        
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
            outHE[HE[HE[he].twin].dest].push_back(he);
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
                double theta;
                Utils::directionAngleInPlane(V[v].pos, V[HE[sortedHE[he]].dest].pos, V[v].n, t1, t2, theta);
                projAngles[he] = theta;
            }
            Utils::doubleListIdxSort(projAngles, sortedHE);
            
            // Rewire based on sorted order
            for (int i = 0; i < sortedHE.size(); i++) {
                int he = sortedHE[i];
                int he_p1 = sortedHE[(i+1)%sortedHE.size()];
                
                int incoming = HE[he_p1].twin;
                HE[he].prev = incoming;
                HE[incoming].next = he;
            }
        }
        return 1;
    }

    // Reset every face based on the halfedges
    int cutmesh::resetFaces() {
        // Throw away all old faces
        F.clear();
        active_f = 0;
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
                seenHE[curr_he] = true;
                HE[curr_he].face = f;
                curr_he = HE[curr_he].next;
                counter++;
            } while (curr_he != he && counter < HE.size());
            if (curr_he != he || heLoop.size() < 3) {
                return -1;
            }
            F.emplace_back();
            active_f++;
            F[f].he = he;
        }
        // Compute areas and normals for all created faces
        computeFNormalsAreas();
        countNumActive();
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
                if ((V[v].projData.elType != 2)) {
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
                active_f--;
            }
        }

        // Loose curves that form no loops.
        // Iterate over halfedges and identify which are without a face and not boundary
        for (int he = 0; he < HE.size(); he++) {
            // NOTE: face index of -1 indicates boundary, meaning that an additional boundary flag set to false is a contradiction
            int twin = HE[he].twin;
            if (HE[he].face == -1 && !HE[he].boundary &&
                HE[twin].face == -1 && !HE[twin].boundary) {
                // Deactivate the associated vertices and edge
                V[HE[he].dest].active = false;
                E[HE[he].edge].active = false;
                HE[he].active = false;
                HE[HE[he].twin].active = false;
            }
        }
        countNumActive();

        return 1;
    }

    // Cuts mesh by "unzipping" along curve vertices
    int cutmesh::cutMesh() {
        // Iterate over each dCNVert
        for (int v = 0; v < V.size(); v++) {
            if (V[v].label != 0 && V[v].active) {  // Only process dCN vertices
                // Gather all adjacent halfedges
                std::vector<int> adjHE = vertAdjHEs(v);
                std::vector<int> HE_SplitIdxs;  // Local HE indexing
                std::vector<int> local_dCNIdxs; // Local HE indexing
                // Find which halfedges to "split" at
                // These include dCN halfedges and boundary edges
                for (int he = 0; he < adjHE.size(); he++) {
                    // Track all partitions that include the boundary
                    if (HE[adjHE[he]].active && ((HE[adjHE[he]].dCN_idx >= 0) || HE[HE[adjHE[he]].twin].boundary)) {
                        HE_SplitIdxs.push_back(he);
                    }
                    // Also track all dCN halfedge Idxs in CCW order
                    if (HE[adjHE[he]].active && HE[adjHE[he]].dCN_idx >= 0) {
                        local_dCNIdxs.push_back(he);
                    }
                }
                if (local_dCNIdxs.empty()) {
                    continue;
                }
                // Find corner index for each
                std::vector<int> cornerIdxs(HE_SplitIdxs.size());   // Global HE indexing
                for (int i = 0; i < HE_SplitIdxs.size(); i++) {
                    int splitLocal = HE_SplitIdxs[i];
                    // Find the corresponding halfedge
                    int next_dCN = local_dCNIdxs[0];
                    for (int idx = 0; idx < local_dCNIdxs.size(); idx++) {
                        int local_i = local_dCNIdxs[idx];
                        if (local_i >= splitLocal) {
                            next_dCN = local_i;
                            break;
                        }
                    }
                    cornerIdxs[i] = HE[adjHE[next_dCN]].twin;
                }

                // For each dCN corner, split the vertex and rewire the halfedge destination vertices accordingly
                std::vector<std::pair<int, int>> newStartBoundaries;
                std::vector<std::pair<int, int>> newEndBoundaries;
                int num_corners = 0;
                for (int i = 0; i < HE_SplitIdxs.size(); i++) {
                    int startHE = adjHE[HE_SplitIdxs[i]];
                    int nextHE = adjHE[HE_SplitIdxs[(i+1)%HE_SplitIdxs.size()]];
                    if (HE[startHE].boundary) { // Skip adding any corners that are boundaries
                        continue;
                    }
                    int v_current;
                    if (num_corners == 0) {
                        v_current = v;
                    } else {
                        v_current = insertVertex(V[v].pos, V[v].n, V[v].label, -1, V[v].projData.elType, V[v].projData.elIdx, V[v].defData.projVector);
                    }
                    // Rewire corner index
                    V[v_current].corner_idx = cornerIdxs[i];
                    V[v_current].he = startHE;

                    // Check if we need to first insert the start boundary
                    std::pair<int, int> startBdy;
                    if (HE[HE[startHE].twin].boundary) {
                        int adj_v = HE[startHE].dest;
                        HE[HE[startHE].twin].dest = v_current;
                        vertPairToHE.erase({adj_v, v});
                        vertPairToHE[{adj_v, v_current}] = HE[startHE].twin;
                        startBdy.first = HE[startHE].twin;
                        startBdy.second = HE[startBdy.first].prev;
                    } else {    // Currently no existing halfedge. Create a new one
                        // Each split that doesn't start at a boundary keeps its last edge and gets a copy for its start edge
                        int new_e = E.size();
                        E.emplace_back();
                        E[new_e].he = startHE;
                        int new_he = HE.size();
                        HE.emplace_back();
                        HE[new_he].twin = startHE;
                        HE[new_he].dest = HE[HE[startHE].twin].dest;
                        HE[new_he].boundary = true;
                        HE[new_he].edge = new_e;
                        HE[startHE].edge = new_e;
                        startBdy.first = new_he;
                        startBdy.second = -1;
                    }
                    newStartBoundaries.push_back(startBdy);
                    int he_current = startHE;
                    int counter = 0;    // safety
                    // Rewire the halfedges to point to the new vert
                    do {
                        int adj_v0 = HE[he_current].dest;
                        vertPairToHE.erase({v, adj_v0});
                        vertPairToHE[{v_current, adj_v0}] = he_current;

                        int he_in = HE[he_current].prev;
                        int adj_v1 = HE[he_in].dest;
                        // Sanity
                        HE[he_in].dest = v_current;
                        if (HE[he_in].dest != v_current) {
                            return -1;
                        }
                        vertPairToHE.erase({adj_v1, v});
                        vertPairToHE[{adj_v1, v_current}] = he_in;
                        
                        he_current = HE[he_in].twin;
                        counter++;
                    } while((he_current != nextHE && !HE[he_current].boundary) && (counter < adjHE.size()));
                    if (he_current != nextHE) {
                        return -1;
                    }

                    // Check if we need to insert the end boundary
                    std::pair<int, int> endBdy; // 
                    if (HE[he_current].boundary) {
                        int adj_v = HE[he_current].dest;
                        vertPairToHE.erase({v, adj_v});
                        vertPairToHE[{v_current, adj_v}] = he_current;
                        endBdy.first = he_current;
                        endBdy.second = HE[he_current].next;
                    } else {    // Need to create a new halfedge
                        int new_he = HE.size();
                        HE.emplace_back();
                        HE[new_he].edge = HE[he_current].edge;
                        HE[new_he].boundary = true;
                        HE[new_he].twin = HE[he_current].twin;
                        HE[new_he].dest = HE[he_current].dest;
                        endBdy.first = new_he;
                        endBdy.second = -1;
                    }
                    newEndBoundaries.push_back(endBdy);
                    num_corners++;
                }

                // Now that all splits are made, rewire by properly rewiring the boundaries
                for (int i = 0; i < newStartBoundaries.size(); i++) {
                    const std::pair<int, int>& currStart = newStartBoundaries[i];
                    const std::pair<int, int>& currEnd = newEndBoundaries[i];

                    HE[currStart.first].next = currEnd.first;
                    HE[currEnd.first].prev = currStart.first;

                    HE[HE[currStart.first].twin].twin = currStart.first;
                    HE[HE[currEnd.first].twin].twin = currEnd.first;

                    // Reconnect the newly formed boundaries
                    // Do both reciprocally for safety
                    if (currStart.second == -1) {
                        const std::pair<int, int>& prevEnd = newEndBoundaries[(i+newEndBoundaries.size()-1)%newEndBoundaries.size()];
                        HE[currStart.first].prev = prevEnd.first;
                        HE[prevEnd.first].next = currStart.first;
                    } else {
                        HE[currStart.first].prev = currStart.second;
                        HE[currStart.second].next = currStart.first;
                    }
                    if (currEnd.second == -1) {
                        const std::pair<int, int>& nextStart = newStartBoundaries[(i+1)%newStartBoundaries.size()];
                        HE[currEnd.first].next = nextStart.first;
                        HE[nextStart.first].prev = currEnd.first;
                    } else {
                        HE[currEnd.first].next = currEnd.second;
                        HE[currEnd.second].prev = currEnd.first;
                    }
                }
            }
        }

        countNumActive();
        return 1;
    }

}   // namespace Mesh