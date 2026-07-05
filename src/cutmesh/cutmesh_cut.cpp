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
            std::vector<Vert> traceList;
            traceList.push_back(V[v0]);
            Eigen::Vector3d direc = (V[v1].pos - V[v0].pos).normalized();
            int success = M->traceGeodesic(V[v0], V[v1], direc, 
                                           V[v0].mesh_elType, V[v0].mesh_elIdx, 
                                           traceList);
            // This is a likely spot for failure, so flag it
            if (success == -1) {
                return -1;
            }
            traceList.push_back(V[v1]);

            // Insert each vert into the cutmesh
            std::vector<int> traceVerts;
            traceVerts.push_back(v0);
            for (int v_idx = 1; v_idx < traceList.size() - 1; v_idx++) {
                traceList[v_idx].label = 2;
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

        // Compute faces
        sortHalfEdges();
        resetFaces();

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

            // Also, this is a good time to label dCN_idx for verts being split!
        }

        // For each outgoing dCN halfedge (above 1), split the vertex as a corner (be sure to add to associated dCNVtoV)
        return 1;
    }

}   // namespace Mesh