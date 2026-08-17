// File for mesh cutting
#include "cutmesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
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
        const std::vector<Polynet::Vert>& dCN_V = dCN->V;
        std::vector<Vert> proj_V(dCN_V.size());
        std::vector<CutData> proj_CD(dCN_V.size());
        // 1. First, create a list of all projections (parallel)
        #pragma omp parallel for
        for (int v = 0; v < dCN_V.size(); v++) {
            // Since we already compute this at the dcurvenet level, grab the data from it
            const Utils::frameData proj = dCN_V[v].proj;
            Eigen::Vector3d n;
            if (proj.proj.elType == 0) {  // Vertex
                n = M->V[proj.proj.elIdx].n;
            } else if (proj.proj.elType == 1) {   // Edge
                n = M->E[proj.proj.elIdx].n;
            } else {    // Face
                n = M->F[proj.proj.elIdx].n;
            }
            proj_V[v] = createVertex(dCN_V[v].pos - proj.offset, n);
            proj_CD[v].label = 1;
            proj_CD[v].corner_idx = -1;
            proj_CD[v].projData.elType = proj.proj.elType;
            proj_CD[v].projData.elIdx = proj.proj.elIdx;
            proj_CD[v].defData.projVector = proj.offset;
        }
        
        // Now, add all vertices into the cutmesh (SEQUENTIAL)
        std::map<int, std::vector<std::pair<double, int>>> edgeMap;
        for (int v = 0; v < proj_V.size(); v++) {
            // First, identify what kind of vertex should be inserted
            int elType = proj_CD[v].projData.elType;
            int elIdx = proj_CD[v].projData.elIdx;
            // If it landed on a vertex, then modify the existing vertex
            if (elType == 0) {
                // If we landed on an existing cut-vertex, then things are bad!
                if (cutData[elIdx].label != 0) {
                    return -1;
                }
                cutData[elIdx].label = 1;
                cutData[elIdx].projData = proj_CD[v].projData;
                cutData[elIdx].defData.projVector = proj_CD[v].defData.projVector;
                dCNVtoV[v] = elIdx;
            } else if (elType == 1) {
                // If it landed on an edge, we need to figure out where exacty to split the edge
                int new_v = insertVertex(proj_V[v], proj_CD[v]);
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
                int new_v = insertVertex(proj_V[v], proj_CD[v]);
                dCNVtoV[v] = new_v;
            }
        }

        // Trace the curves (SEQUENTIAL, due to geodesics)
        const std::vector<Polynet::Edge>& dCN_E = dCN->E;
        const std::vector<Polynet::HalfEdge>& dCN_HE = dCN->HE;
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
            std::vector<Utils::projData> traceProjData;
            traceVerts.push_back(V[v0]);
            traceProjData.push_back(cutData[v0].projData);
            int depth = 0;
            Eigen::Vector3d direc = (V[v1].pos - V[v0].pos).normalized();
            //Eigen::Vector3d direc = (dCN_V[dCN_v1].pos - dCN_V[dCN_v0].pos).normalized();
            std::cout << "Starting trace from vert " << v0 << " (proj type: " << cutData[v0].projData.elType << ") to " << v1 << " (proj type: " << cutData[v1].projData.elType << ")" << std::endl;
            /*
            std::cout << "Initial direction: " << direc[0] << ", "
                                        << direc[1] << ", "
                                        << direc[2] << std::endl;
            */
            int success = M->traceGeodesic(V[v0], cutData[v0].projData, V[v1], cutData[v1].projData, direc,
                                           cutData[v0].projData,
                                           traceVerts, traceProjData, depth);
            // This is a likely spot for failure, so flag it
            if (success != 1) {
                std::cout << "Trace failed. Trying opposite direction." << std::endl;
                // Check the opposite direction to catch initial directional error
                // TODO: How necessary is this?
                depth = 0;
                traceVerts.clear();
                traceProjData.clear();
                traceVerts.push_back(V[v0]);
                traceProjData.push_back(cutData[v0].projData);
                int success_opposite = M->traceGeodesic(V[v0], cutData[v0].projData, V[v1], cutData[v1].projData, -1 * direc,
                                           cutData[v0].projData,
                                           traceVerts, traceProjData, depth);
                if (success_opposite != 1) {
                    std::cout << "Trace failed again. Quitting." << std::endl;
                    return -1;
                }
            }
            traceVerts.push_back(V[v1]);
            traceProjData.push_back(cutData[v1].projData);

            // Insert each vert into the cutmesh
            std::vector<int> traceList;
            traceList.push_back(v0);
            for (int v_idx = 1; v_idx < traceVerts.size() - 1; v_idx++) {
                CutData traceCD;
                traceCD.label = 2;
                traceCD.projData = traceProjData[v_idx];
                // TODO: Is this a safe thing to do?
                // i.e., the projVector taken using the midpoint of the associated halfedge, thus being deterministic without
                // Requiring us to subdivide the dCN edge OR change the projection estimate formula (sort of)
                traceCD.defData.projVector = ((dCN_V[dCN_v1].pos + dCN_V[dCN_v0].pos) / 2) - traceVerts[v_idx].pos;
                // First, compute an estimated curvenet position so we can take the difference
                if (traceProjData[v_idx].elType == 0) {  // Check if we are on a vertex
                    if (cutData[traceProjData[v_idx].elIdx].label != 0) { // If we hit a vertex that is already assigned, then quit
                        return -1;
                    }
                    cutData[traceProjData[v_idx].elIdx].label = 2;
                    cutData[traceProjData[v_idx].elIdx].projData = traceCD.projData;
                    cutData[traceProjData[v_idx].elIdx].defData = traceCD.defData;
                    traceList.push_back(traceProjData[v_idx].elIdx);
                } else if (traceProjData[v_idx].elType == 1) { // We must be on an edge
                    int new_v = insertVertex(traceVerts[v_idx], traceCD);
                    // Check if the edge was already split. If so, find where to split it.
                    int orig_e = cutData[new_v].projData.elIdx;   // source mesh edge index
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
                    V[v_start].he = he0;
                    V[v_next].he = he1;
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
        deactivateIsolatedCuts();
        std::cout << "Isolated cuts deactivated." << std::endl;
        
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
            V[v].he = sortedHE[0];
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
                if (HE[HE[he].twin].dest != v) {
                    std::cout << "Mesh halfedges violate twin rule." << std::endl;
                    return -1;
                }
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
                if ((cutData[v].projData.elType != 2)) {
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
        int max_v = V.size();
        // Iterate over each dCNVert
        for (int v = 0; v < max_v; v++) {
            if (cutData[v].label == 0 || !V[v].active) { // Only process dCN vertices
                continue;
            }
            // Gather all adjacent halfedges
            std::vector<int> adjHE = vertAdjHEs(v);
            std::vector<int> adjV = vertAdjVerts(v);
            std::vector<int> incomingHE(adjHE.size());
            std::vector<int> HE_SplitIdxs;  // Local HE indexing
            std::vector<int> local_dCNIdxs; // Local HE indexing
            for (int he = 0; he < adjHE.size(); he++) { // Fill in the incoming halfedges
                incomingHE[he] = HE[adjHE[he]].twin;
            }
            // Find which halfedges to "split" at
            // These include dCN halfedges and boundary edges
            for (int he = 0; he < adjHE.size(); he++) {
                // Track all partitions that are embedded and NOT boundary, or the starts of boundaries
                if (HE[adjHE[he]].active && (((HE[adjHE[he]].dCN_idx >= 0) && !HE[adjHE[he]].boundary) || 
                    HE[HE[adjHE[he]].twin].boundary)) {
                    HE_SplitIdxs.push_back(he);
                }
                // Also track all dCN halfedge Idxs in CCW order
                if (HE[incomingHE[he]].active && HE[incomingHE[he]].dCN_idx >= 0) {
                    local_dCNIdxs.push_back(he);
                }
            }
            if (local_dCNIdxs.empty()) {
                continue;
            }
            // Find corner index for each
            std::vector<int> cornerIdxs(HE_SplitIdxs.size());   // Global HE indexing
            for (int i = 0; i < HE_SplitIdxs.size(); i++) {     // Iterate over each new corner
                int splitLocal = HE_SplitIdxs[i];
                // Find the corresponding halfedge
                int next_dCN = local_dCNIdxs[0];
                for (int idx = 0; idx < local_dCNIdxs.size(); idx++) {  // Find this corner's corner_idx
                    int local_i = local_dCNIdxs[idx];
                    if (local_i >= splitLocal) {
                        next_dCN = local_i;
                        break;
                    }
                }
                cornerIdxs[i] = incomingHE[next_dCN];
            }

            // For each dCN corner, split the vertex and rewire the halfedge destination vertices accordingly
            std::vector<std::pair<int, int>> newStartBoundaries;
            std::vector<std::pair<int, int>> newEndBoundaries;
            for (int i = 0; i < HE_SplitIdxs.size(); i++) {
                int startHE = adjHE[HE_SplitIdxs[i]];
                int nextHE = adjHE[HE_SplitIdxs[(i+1)%HE_SplitIdxs.size()]];
                int startLocal = HE_SplitIdxs[i];
                int endLocal = HE_SplitIdxs[(i + 1) % HE_SplitIdxs.size()];

                int num_outgoing = endLocal - startLocal;
                if (num_outgoing <= 0) {
                    num_outgoing += adjHE.size();
                }
                int v_current;
                if (i == 0) {
                    v_current = v;
                } else {
                    CutData splitCD;
                    splitCD.label = cutData[v].label;
                    splitCD.projData = cutData[v].projData;
                    splitCD.defData.projVector = cutData[v].defData.projVector;
                    v_current = insertVertex(V[v].pos, V[v].n, splitCD);
                }
                // Rewire corner index
                cutData[v_current].corner_idx = cornerIdxs[i];
                V[v_current].he = startHE;

                // Check if we need to first insert the start boundary
                std::pair<int, int> startBdy = {-1, -1};
                if (HE[HE[startHE].twin].boundary) {
                    int adj_v = HE[startHE].dest;
                    HE[HE[startHE].twin].dest = v_current;
                    vertPairToHE.erase({adj_v, v});
                    vertPairToHE[{adj_v, v_current}] = HE[startHE].twin;
                    startBdy.first = HE[startHE].twin;
                    startBdy.second = HE[startBdy.first].prev;
                } else {    // Currently no existing halfedge. Create a new one
                    // Each split that doesn't start at a boundary keeps its last edge and gets a copy for its start edge
                    int adj_v = HE[startHE].dest;
                    int new_e = E.size();
                    E.emplace_back();
                    E[new_e].he = startHE;
                    int new_he = HE.size();
                    HE.emplace_back();
                    HE[new_he].twin = startHE;
                    HE[new_he].dest = v_current;
                    HE[new_he].boundary = true;
                    HE[new_he].edge = new_e;
                    HE[startHE].edge = new_e;
                    vertPairToHE.erase({adj_v, v});
                    vertPairToHE[{adj_v, v_current}] = new_he;
                    startBdy.first = new_he;
                    startBdy.second = -1;
                }
                newStartBoundaries.push_back(startBdy);
                int counter = 0;    // safety
                int local_final_HE = endLocal;
                // Rewire the halfedges to point to the new vert
                for (int j = 0; j < num_outgoing; j++) {
                    int local = (startLocal + j) % adjHE.size();
                    int nextLocal = (local + 1) % adjHE.size();

                    int he_current = adjHE[local];
                    int adj_v0 = HE[he_current].dest;

                    // he_current was v -> adj_v0.
                    // It should become v_current -> adj_v0.
                    vertPairToHE.erase({v, adj_v0});
                    vertPairToHE[{v_current, adj_v0}] = he_current;

                    // If this outgoing halfedge is itself a boundary, this sector ends here.
                    // Do not try to use the next incoming twin.
                    if (HE[he_current].boundary) {
                        local_final_HE = local;
                        break;
                    }

                    // In the sorted ring convention, the incoming halfedge closing this wedge
                    // is the twin of the next outgoing halfedge.
                    int he_prev = HE[adjHE[nextLocal]].twin;

                    // he_prev should originally be adj_v1 -> v.
                    int adj_v1 = HE[adjHE[nextLocal]].dest;

                    if (HE[he_prev].dest != v) {
                        std::cout << "Reached corner " << i + 1 << " of " << HE_SplitIdxs.size()
                                << "; Expected incoming halfedge to end at split vertex." << std::endl;

                        std::cout << "  v = " << v
                                << ", v_current = " << v_current << std::endl;
                        std::cout << "  local = " << local
                                << ", nextLocal = " << nextLocal
                                << ", startLocal = " << startLocal
                                << ", endLocal = " << endLocal
                                << ", num_outgoing = " << num_outgoing << std::endl;
                        std::cout << "  he_current = " << he_current
                                << ", he_prev = " << he_prev
                                << ", adjHE[nextLocal] = " << adjHE[nextLocal] << std::endl;
                        std::cout << "  HE[he_prev].dest = " << HE[he_prev].dest
                                << ", expected = " << v << std::endl;

                        return -1;
                    }

                    // he_prev was adj_v1 -> v.
                    // It should become adj_v1 -> v_current.
                    HE[he_prev].dest = v_current;

                    vertPairToHE.erase({adj_v1, v});
                    vertPairToHE[{adj_v1, v_current}] = he_prev;
                }

                // Check if we need to insert the end boundary
                int final_outgoing = adjHE[local_final_HE];
                std::pair<int, int> endBdy; // 
                if (HE[final_outgoing].boundary) {
                    int adj_v = adjV[local_final_HE];
                    vertPairToHE.erase({v, adj_v});
                    vertPairToHE[{v_current, adj_v}] = final_outgoing;
                    endBdy.first = final_outgoing;
                    endBdy.second = HE[final_outgoing].next;
                } else {    // Need to create a new halfedge
                    int adj_v = adjV[local_final_HE];
                    int new_he = HE.size();
                    HE.emplace_back();
                    HE[new_he].edge = HE[final_outgoing].edge;
                    HE[new_he].boundary = true;
                    HE[new_he].twin = HE[final_outgoing].twin;
                    HE[new_he].dest = adj_v;
                    vertPairToHE.erase({adj_v, v});

                    vertPairToHE[{v_current, adj_v}] = new_he;
                    endBdy.first = new_he;
                    endBdy.second = -1;
                }
                newEndBoundaries.push_back(endBdy);
            }

            // Now that all splits are made, rewire by properly rewiring the boundaries
            for (int i = 0; i < newStartBoundaries.size(); i++) {
                const std::pair<int, int>& currStart = newStartBoundaries[i];
                const std::pair<int, int>& currEnd = newEndBoundaries[i];
                const std::pair<int, int>& nextStart = newStartBoundaries[(i+1)%newStartBoundaries.size()];
                const std::pair<int, int>& prevEnd = newEndBoundaries[(i+newEndBoundaries.size()-1)%newEndBoundaries.size()];
                // Connect the start boundary to the end boundary
                HE[currStart.first].next = currEnd.first;
                HE[currEnd.first].prev = currStart.first;
                // Make sure twins are set correctly
                HE[HE[currStart.first].twin].twin = currStart.first;
                HE[HE[currEnd.first].twin].twin = currEnd.first;

                // Reconnect the newly formed boundaries
                // Do both reciprocally for safety
                if (currStart.second == -1 && prevEnd.second == -1) {   // created a new boundary
                    HE[currStart.first].prev = prevEnd.first;
                    HE[prevEnd.first].next = currStart.first;
                } else if (currStart.second >= 0 && prevEnd.second >= 0) {  // Joining up existing boundary halfedges
                    HE[currStart.first].prev = currStart.second;
                    HE[currStart.second].next = currStart.first;
                } else {        // Mismatch!
                    std::cout << "Mismatched adjacent start/prev boundaries" << std::endl;
                    return -1;
                }
                if (currEnd.second == -1 && nextStart.second == -1) { // Using an existing boundary
                    HE[currEnd.first].next = nextStart.first;
                    HE[nextStart.first].prev = currEnd.first;
                } else if (currEnd.second >= 0 && nextStart.second >= 0) {  // Joining up existing boundary vertices
                    HE[currEnd.first].next = currEnd.second;
                    HE[currEnd.second].prev = currEnd.first;
                } else {    // Mismatch!
                    std::cout << "Mismatched adjacent end/next boundaries" << std::endl;
                    return -1;
                }
            }
        }

        countNumActive();
        return 1;
    }

}   // namespace Mesh