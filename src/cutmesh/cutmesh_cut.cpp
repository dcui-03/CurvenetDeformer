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
        dCNVtoV.clear();
        dCNHEtoHE.clear();
        
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
            proj_V[v] = createVertex(proj, n, 1, v, elType, elIdx, pos-proj);
        }
        
        // Now, add all vert and face vertices into the cutmesh (SEQUENTIAL)
        std::map<int, std::vector<Vert>> edgeMap;
        for (int v = 0; v < proj_V.size(); v++) {
            // First, identify what kind of vertex should be inserted
            int elType = proj_V[v].mesh_elType;
            int elIdx = proj_V[v].mesh_elIdx;
            // If it landed on a vertex, then modify the existing vertex
            if (elType == 0) {
                V[elIdx] = proj_V[v];
                dCNVtoV[v].push_back(elIdx);
            } else if (elType == 1) {

                // If it landed on an edge, then create the vertex and split the edge
                int new_v = insertVertex(proj_V[v]);
                splitEdge(elIdx, new_v);
                dCNVtoV[v].push_back(new_v);
            } else {
                // If it landed on a face, then simply insert the vertex
                int new_v = insertVertex(proj_V[v]);
                dCNVtoV[v].push_back(new_v);
            }
        }

        // Add in all edge splits
        for (const std::pair<int, std::vector<Vert>> split : edgeMap) {
            int e = split.first;
            std::vector<Vert> splitVerts = split.second;
            // Sort by t-value
            // TODO
            // Split edges sequentially by t-value
        }

        // Trace the curves (SEQUENTIAL, due to geodesics)
        const std::vector<DCurvenet::Edge>& dCN_E = dCN->E;
        const std::vector<DCurvenet::HalfEdge>& dCN_HE = dCN->HE;
        for (int e = 0; e < dCN_E.size(); e++) {
            // For this edge, connect the two associated vertices in the cutmesh
            int dCN_he0 = dCN_E[e].he;
            int dCN_he1 = dCN_HE[dCN_he0].twin;
            int v0 = dCNVtoV[dCN_HE[dCN_he1].dest][0];
            int v1 = dCNVtoV[dCN_HE[dCN_he0].dest][0];

            // If edge already exists, modify the existing edge
            if (vertPairToHE.find(std::make_pair(v0, v1)) != vertPairToHE.end()) {
                int he0 = vertPairToHE[std::make_pair(v0, v1)];
                int he1 = HE[he0].twin;
                HE[he0].dCN_idx = dCN_he0;
                HE[he1].dCN_idx = dCN_he1;
                
                dCNHEtoHE[dCN_he0].push_back(he0);
                dCNHEtoHE[dCN_he0].push_back(he1);
                continue;
            }
            // Classify based on properties
            if (V[v0].mesh_elType == 2 && V[v1].mesh_elType == 2) {  // Case 1: Two face verts
                // Subcase 1: Both are on same face (just draw the edge)
                if (V[v0].mesh_elIdx == V[v1].mesh_elIdx) {
                    insertEdge(v0, v1, dCN_he0, dCN_he1);
                } else { // Subcase 2: Not on same face --> Trace Geodesic
                    std::vector<Vert> traceList;
                    traceList.push_back(V[v0]);
                    Eigen::Vector3d direc = (V[v1].pos - V[v0].pos).normalized();
                    int success = M->traceGeodesic(V[v0], V[v1], direc, traceList, dCN_he0, dCN_he1, 2);
                    // This is a likely spot for failure, so flag it
                    if (success == -1) {
                        return -1;
                    }
                    traceList.push_back(V[v1]);

                    // Insert each vert into the cutmesh
                    std::vector<int> traceVerts;
                    traceVerts.push_back(v0);
                    for (int v_idx = 1; v_idx < traceList.size() - 1; v_idx++) {
                        if (traceList[v_idx].mesh_elType == 0) {  // Check if we are on a vertex
                            V[traceList[v_idx].mesh_elIdx] = traceList[v_idx];
                            traceVerts.push_back(traceList[v_idx].mesh_elIdx);
                        } else { // Check if we are on an edge
                            int new_v = insertVertex(traceList[v_idx]);
                            // TODO: Check if the edge was already split. If so, find where to split it.
                            traceVerts.push_back(new_v);
                        }
                    }
                    traceVerts.push_back(v1);

                    // Connect the inserted vertices
                    for (int v_idx = 1; v_idx < traceList.size(); v_idx++) {
                        int e = insertEdge(traceVerts[v_idx - 1], traceVerts[v_idx], dCN_he0, dCN_he1);
                        dCNHEtoHE[dCN_he0].push_back(E[e].he);
                        dCNHEtoHE[dCN_he1].push_back(HE[E[e].he].twin);
                    }
                }
            } else {    // Case 2: Everything else (verts already inserted)
                // Add shared edge
                int e = insertEdge(v0, v1, dCN_he0, dCN_he1);
                dCNHEtoHE[dCN_he0].push_back(E[e].he);
                dCNHEtoHE[dCN_he0].push_back(HE[E[e].he].twin);
            }
        }

        // Compute faces
        sortHalfEdges();
        resetFaces();

        return 1;
    }

    // Properly orient halfedges based on local tangent planes
    int cutmesh::sortHalfEdges() {
        // Determine incoming halfedges for all vertices
        std::vector<std::vector<int>> incomingVHalfEdges(V.size());
        for (int he = 0; he < HE.size(); he++) {
            if (!HE[he].active) {
                continue;
            }
            incomingVHalfEdges[HE[HE[he].dest].twin].push_back(he);
        }

        // Iterate over vertices and sort their halfedges based on projection to the tangent plane
        for (int v = 0; v < V.size(); v++) {
            if (!V[v].active) {
                continue;
            }
            std::vector<Eigen::Vector3d> adjVecs(incomingVHalfEdges[v].size());
            for (int he = 0; he < incomingVHalfEdges[v].size(); he++) {
                adjVecs[he] = (V[HE[he].dest].pos - V[v].pos).normalized();
            }
            // Project all onto local tangent plane and sort
            // TODO
            
            // Rewire based on sorted order
            std::vector<int> sortedHalfEdges(adjVecs.size());
            for (int i = 0; i < sortedHalfEdges.size(); i++) {
                int he = sortedHalfEdges[i];
                int he_p1 = sortedHalfEdges[(i+1)%sortedHalfEdges.size()];
                
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
                HE[curr_he] = true;
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

    // Cuts mesh by "unzipping" along curve vertices
    int cutmesh::cutMesh() {
        // Iterate over each dCNVert
        for (int v = 0; v < V.size(); v++) {
            if (V[v].label == 0) {
                continue;
            }
            std::vector<int> adjHEs = vertAdjHEs(v);
            std::vector<bool> is_CNHE(adjHEs.size());
            int num_split = 0;
            // Check which halfedges are from the curvenet
            for (int he : adjHEs) {
                if (HE[adjHEs[he]].dCN_idx != -1) {
                    is_CNHE[he] = true;
                    num_split++;
                }
            }
            // For each outgoing dCN HE whose twin is also a dCN HE, create a duplicate vertex,
            // BUT only do this if the number of outgoing dCN HE's is 
            // TODO
            // To "split" the halfedges, simply reset the dest vertex of the incoming halfedge and reset that vertex's pointed-to halfedge
            // BE CAREFUL AT BOUNDARIES! Do NOT split if a dCN edge is on a boundary
        }

        // For each outgoing dCN halfedge (above 1), split the vertex as a corner (be sure to add to associated dCNVtoV)
        return 1;
    }

    // Recursively traces a geodesic from an initial position to a "goal" point
    // NOTE: dCN_idx here is a halfedge index of the discrete curve network
    // TODO: Do we want to use the vector method? This may help with straightest geodesics, but may also become unreliable
    //       However, the opposite is true about projection at each step... we are guaranteed to reach the goal, but may end up not tracing 
    //       an accurate geodesic... (accuracy is not a big problem if sufficient curve samples are taken).
    //       But note also that projection at each step also prevents drift due to vertex-snapping and non-planar faces
    //       A compromise is to compute the next vector direction at the end of each call, then pass it in to the next call, which lets you pre-compute an accurate very first direction.
    //       Maybe also store a nextType and nextIdx as inputs? These might help speed things up
    
    int cutmesh::traceGeodesic(int start, int end, int dCN_idx0, int dCN_idx1, Eigen::Vector3d direc, bool reDirec, bool fast) {
        // 1. Identify start adjacent faces and end adjacent faces
        //    If the start or end vertex has no adjacent faces (i.e., halfedges) then it must be new and in a face. Take its input idx instead.
        //    NOTE: vertAdjFaces() already returns only unique adjacent faces
        double eps = 1e-3;
        std::vector<int> start_adjF = vertAdjFaces(start);
        std::vector<int> end_adjF = vertAdjFaces(end);
        std::vector<int> sharedFaces;
        for (int i = 0; i < start_adjF.size(); i++) {
            for (int j = 0; j < end_adjF.size(); j++) {
                if ((start_adjF[i] != -1) && (start_adjF[i] == end_adjF[j])) {  // Shared face is not boundary
                    sharedFaces.push_back(start_adjF[i]);
                }
            }
        }
        // Case 1: Both vertices are isolated/new
        // 2. If both vertices have no adjacent faces, then they must be newly initialized with no
        //    adjacent edges. Check the passed-in starting and ending faces to see if we can go ahead and connect them
        if (start_adjF.size() == 0 && end_adjF.size() == 0) {
            if ((startF == -1) && (endF == -1)) {   // Unlikely, unless func is called incorrectly
                return -1;
            } else if (startF == -1) {  // Also unlikely, but simply draw an edge between them
                int new_e = insertEdge(endF, start, end, dCN_idx0, dCN_idx1, true);
                return 1;
            } else if (endF == -1) {
                int new_e = insertEdge(startF, start, end, dCN_idx0, dCN_idx1, true);
                return 1;
            } else if (startF == endF) {   // New isolated face-points on the same face. Insert the edge
                int new_e = insertEdge(startF, start, end, dCN_idx0, dCN_idx1, true);
                return 1;
            } else {
                // TODO: Special case. Project each point onto each of their corresp. face planes and compute next geodesic intersection

            }
        }
        // Case 2: 2 faces are shared by the vertex pair. Note for a manifold surface, 2 should be the maximum
        else if (sharedFaces.size() >= 2) { // If both the start and end are on vertices and are neighbors, then simply modify the existing edge.
            int he0 = vertPairToHE[std::make_pair(start, end)];
            int he1 = HE[he0].twin;
            HE[he0].dCN_idx = dCN_idx0;
            HE[he1].dCN_idx = dCN_idx1;
            return 1;   // success
        } 
        // Case 3: If they share an adjacent face or are visible on that face, then draw an edge on that face and return
        else if (sharedFaces.size() == 1) {
            if (fast) { // FAST VERSION: Simple shared face check (better is all faces are convex/near convex)
                int new_e = insertEdge(sharedFaces[0], start, end, dCN_idx0, dCN_idx1, true);
                return 1;
            }
        }
        
        // If none of these tests are passed, then defer to to tracing the geodesic
        // 3. If they do not pass the neighbor test, then grab the vector from the start to the end and project it onto the tangent plane of ANY face adjacent to the start
        if (reDirec) {  // Recompute global walk direction
            direc = (V[end].pos - V[start].pos).normalized();
        }

        // 4. Test if the vector is coincident to any adjacent edge (within a tiny tolerance). If so, snap to that edge
        bool next_reDirec = true;
        int curr_face = -1;  // Face that we are walking on next
        std::vector<int> start_adjHE = vertAdjHEs(start);
        for (int he0_idx = 0; he0_idx < start_adjHE.size(); he0_idx++) {
            int he0 = start_adjHE[he0_idx];
            int he1 = start_adjHE[(he0_idx + 1) % start_adjHE.size()];
            Eigen::Vector3d he0_vec = (V[HE[he0].dest].pos - V[start].pos).normalized();
            Eigen::Vector3d he1_vec = (V[HE[he1].dest].pos - V[start].pos).normalized();
            int f = HE[he0].face;
            if (f == -1) {  // Special handling of boundaries
                // TODO: We need to check if our direction lies between the two faces and snap to a boundary edge in this case
                continue;
            }
            // Project direc onto each adj face
            Eigen::Vector3d proj_direc;
            double valid_proj = Utils::projectVectorOntoTangentPlane(F[f].n, direc, proj_direc);
            if (valid_proj == -1.0) {
                continue;
            }
            // If we are nearly coincident with an edge, just modify that edge
            if (proj_direc.dot(he0_vec) <= 1-eps) {
                int he0_twin = HE[he0].twin;
                HE[he0].dCN_idx = dCN_idx0;
                HE[he0_twin].dCN_idx = dCN_idx1;
                int new_start = HE[he0].dest;
                return traceGeodesic(new_start, f, end, endF, dCN_idx0, dCN_idx1, he0_vec, true, fast);
            } else if (proj_direc.dot(he1_vec) >= 1-eps) {  // Slightly redundant, but safer for boundary vertices
                int he1_twin = HE[he1].twin;
                HE[he1].dCN_idx = dCN_idx0;
                HE[he1_twin].dCN_idx = dCN_idx1;
                int new_start = HE[he1].dest;
                return traceGeodesic(new_start, f, end, endF, dCN_idx0, dCN_idx1, he1_vec, true, fast);
            }
        }   
        // Failed to find a valid next direction
        if (curr_face == -1) {
            return -1;
        }
        // Project onto each adjacent face to get the "from" direction/face, then use 
        // . Otherwise if the start is a vertex/edge, then check to make sure that we are on the right adjacent face by testing the cross product between the halfedge and the vector
        //    If the cross product is positive when dotted with the face normal, then we are pointing away from the face and should choose a different one.
        //    If the cross product is negative when dotted with the face normal, then we are pointing into the face and can continue
        // 7. Then, project onto the face's Newell plane (or check height) and raycast until a non-starting point intersection is reached.
        // 8. If the intersection is close to a vertex, snap to that vertex and apply modifications to that vertex
        // 9. If landed on an edge, add a new vertex by splitting the edge
        // 10. Add an edge between the start vertex and this new vertex
        // 11. For the new halfedges produced by this edge, add them to the corresponding dCN Halfedge to mesh halfedge map
        // 11. Call recursively, with the current and new vertices as the new start and end
    }

    // Checks whether two points on the mesh are "visible"
    // Includes a fast version (simple mesh object adjacency check) and a slow version (true "visibility").

}   // namespace Mesh