#include "cutmesh.hpp"

#include "dcurvenet/pdcurvenet.hpp"
#include "dcurvenet/components/dvert.hpp"
#include "utils/utils.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <map>
#include <vector>
#include <utility>


namespace Mesh {

    bool cutmesh::adjacencyCheck(const DCurvenet::dvert& start, const DCurvenet::dvert& target, std::vector<int>& sharedFaceList) {
        std::vector<int> adjFacesStart;
        std::vector<int> adjFacesTarget;
        // 1. Grab the adjacent faces of the start vertex
        std::pair<int, int> startProj = start.getProjection();
        if (startProj.first == 2) {   // Vert
            heMesh.vertex_face_neighbors(startProj.second, adjFacesStart);
        } else if (startProj.first == 1) {  // Edge
            heMesh.edgeFaces(startProj.second, adjFacesStart);
        } else {    // Face
            adjFacesStart.push_back(startProj.second);
        }

        // 2. Grab the adjacent faces of the end vertex
        std::pair<int, int> targetProj = target.getProjection();
        if (targetProj.first == 2) {   // Vert
            heMesh.vertex_face_neighbors(targetProj.second, adjFacesTarget);
        } else if (targetProj.first == 1) {  // Edge
            heMesh.edgeFaces(targetProj.second, adjFacesTarget);
        } else {    // Face
            adjFacesTarget.push_back(targetProj.second);
        }

        // Check shared faces
        for (int fi = 0; fi < adjFacesStart.size(); fi++) {
            for (int fj = 0; fj < adjFacesTarget.size(); fj++) {
                if ((adjFacesStart[fi] == adjFacesTarget[fj]) && (adjFacesStart[fi] != -1)) {   // Screen out boundary edges
                    sharedFaceList.push_back(adjFacesStart[fi]);
                }
            }
        }
        // If no shared faces, then return false
        if (sharedFaceList.size() == 0) {
            return false;
        }

        // Special exceptions
        // VV Case: both lie on vertices and share >0 faces. Return true (i.e., both verts are connected by an edge)
        if ((startProj.first == 2) && (targetProj.first == 2) &&
            heMesh.edgeIdxFromVerts(startProj.second, targetProj.second) != -1) {
            return true;
        }

        // FF Case: both lie on the same face AND face is convex. Return true
        if ((startProj.first == 0) && (targetProj.first == 0) && Convex[startProj.second]) {
            return true;
        }

        // VE Case: Check if edge contains the vertex
        if ((startProj.first == 1) && (targetProj.first == 0)) {
            std::vector<int> eVerts;
            heMesh.edgeVertices(startProj.second, eVerts);
            if (eVerts[0] == targetProj.second || eVerts[1] == targetProj.second) {
                return true;
            }
        }
        // EV case
        if ((startProj.first == 0) && (targetProj.first == 1)) {
            std::vector<int> eVerts;
            heMesh.edgeVertices(targetProj.second, eVerts);
            if (eVerts[0] == startProj.second || eVerts[1] == startProj.second) {
                return true;
            }
        }

        // EE case
        if ((startProj.first == 1) && (targetProj.first == 1) &&
            (startProj.second == targetProj.second)) {
            return true;
        }

        // If didn't pass any special exceptions, then we can't be certain and need to run more tests
        return false;
    }

    // Recursively compute straightest geodesic between any two points
    // start is the current start point, target is the goal point,
    // Direc is the movement direction, where el_type and el_idx is what element we should walk on
    bool cutmesh::computeStraightestGeodesic(DCurvenet::dvert& start, DCurvenet::dvert& target,
                                    Eigen::Vector3d direc, int el_type, int el_idx,
                                    std::vector<DCurvenet::dvert>& dVertList, int& depth) {
        // Recursion condition: Check if we have reached the end element and return true
        // 1. Check if we are already adjacent to the target
        std::vector<int> sharedFaces;
        if (adjacencyCheck(start, target, sharedFaces)) {    // We can reach the target using the current start
            return true;
        }

        // 2. If we can't immediately rule them out, then walk along any suspect faces
        DCurvenet::dvert intersection; // Create a new dvert object to store the intersection data
        bool recalibrate = false;
        int associatedHalfEdge = -1;    // The most relevant edge (ex. the intersected edge, the edge walked along, the CCW edge to the vertex intx)
        // First test the shared faces (if any)
        if (sharedFaces.size() > 0) {
            double nearest_t = std::numeric_limits<double>::infinity();
            for (int f = 0; f < sharedFaces.size(); f++) {
                // Project face onto Newell plane centered at start location
                std::vector<Eigen::Vector3d> fVerts(F[sharedFaces[f]].size());
                for (int v = 0; v < F[el_idx].size(); v++) {
                    fVerts[v] = V[F[sharedFaces[f]][v]];
                }
                double ray_t, u, theta;
                int local_idx, intx_type;
                bool hit = Utils::computeFaceIntersectionTarget(fVerts, fNormals[sharedFaces[f]], start.getPosition(), target.getPosition(), 
                                        ray_t, intx_type, local_idx, u, theta, meanE*1e-4);
                if (hit) {  // If we hit the target first, then we can terminate immediately
                    return true;
                }
                // Otherwise, we need to decide which result to keep and initialize the intersection with that result
                if (ray_t < nearest_t) {
                    Eigen::Vector3d tempNormal;
                    int intersection_idx;
                    if (intx_type == 1) {
                        recalibrate = false;
                        // Get edge from endpoints
                        intersection_idx = heMesh.edgeIdxFromVerts(F[sharedFaces[f]][local_idx],  F[sharedFaces[f]][(local_idx + 1)%fVerts.size()]);
                        tempNormal = eNormals[intersection_idx];
                    } else {
                        recalibrate = true;
                        intersection_idx = F[sharedFaces[f]][local_idx];
                        tempNormal = vNormals[intersection_idx];
                    }
                    intersection.setPosition(fVerts[local_idx] + u*(fVerts[(local_idx + 1) %fVerts.size()] - fVerts[local_idx]));
                    intersection.setProjection(intx_type, F[sharedFaces[f]][local_idx]);
                    // intersection.setNormal(tempNormal, true);    // No need to set normals for these temporary dverts
                    associatedHalfEdge = heMesh.directed_edge2he_index(F[sharedFaces[f]][local_idx],  F[sharedFaces[f]][(local_idx + 1)%fVerts.size()]);
                }
            }
        } else { // Otherwise, the target is NOT visible, so we need to walk along the specified directions
            if (el_type == 0) {    // i.e., walk along specified face
                // Get current face and compute intersection data
                std::vector<Eigen::Vector3d> fVerts(F[el_idx].size());
                for (int v = 0; v < F[el_idx].size(); v++) {
                    fVerts[v] = V[F[el_idx][v]];
                }
                double ray_t, u, theta;
                int local_idx, intx_type;
                Utils::computeFaceIntersection(fVerts, fNormals[el_idx], start.getPosition(), direc, 
                                        ray_t, intx_type, local_idx, u, theta, meanE*1e-4);
                Eigen::Vector3d tempNormal;
                int intersection_idx;
                if (intx_type == 1) {   // Landed on an edge
                    recalibrate = false;
                    // Get edge from endpoints
                    intersection_idx = heMesh.edgeIdxFromVerts(F[el_idx][local_idx],  F[el_idx][(local_idx + 1)%fVerts.size()]);
                    tempNormal = eNormals[intersection_idx];
                } else {    // Landed on a vertex
                    recalibrate = true;
                    intersection_idx = F[el_idx][local_idx];
                    tempNormal = vNormals[intersection_idx];
                }
                intersection.setPosition(fVerts[local_idx] + u*(fVerts[(local_idx + 1) %fVerts.size()] - fVerts[local_idx]));
                intersection.setProjection(intx_type, intersection_idx);
            } else {        // i.e., walking along edge
                recalibrate = true;
                std::vector<int> adjVerts;
                heMesh.edgeVertices(el_idx, adjVerts);
                // Grab the endpoint as the next start, but check which way our walk is pointing
                if (direc.dot(V[adjVerts[1]] - V[adjVerts[0]]) >= 0.0) {
                    intersection.setPosition(V[adjVerts[1]]);
                    intersection.setProjection(2, adjVerts[1]);
                } else {
                    intersection.setPosition(V[adjVerts[0]]);
                    intersection.setProjection(2, adjVerts[0]);
                }
            }
        }
        // Append intersection point to the list
        dVertList.push_back(intersection);

        // Compute new starting direction and element
        std::pair<int, int> intx_data = intersection.getProjection();
        int new_el_type, new_el_idx;
        if (intx_data.first == 2) {   // Landed on a vertex
            if (recalibrate) {  // Recompute direction based on target
                // Recalibrate by projecting (intersection - end) onto intersection element's tangent plane
                // If direction is almost identical to the current normal, then defer to standard algorithm
                bool valid = true;

                // TODO: Also check if our new direction ends up point toward a boundary face. If so, terminate.
                if (-1) {
                    return false;
                }
                // If valid direction, recurse using the recalibrated direction
                if (valid) {
                    // recurse
                    return true;
                }
                // Otherwise, fall back on the standard algorithm
            }

            // Standard algorithm: Get the CCW fan of the current vertex (including boundary face)
            std::vector<int> neighbors;
            heMesh.vertex_vertex_neighbors(intx_data.second, neighbors);
            // Get total Gaussian curvature using fan. (INCLUDING BOUNDARY)
            double k = 0.0;
            for (int n = 0; n < neighbors.size(); n++) {
                k += Utils::vectorAngle(V[intx_data.second], V[neighbors[n]],
                                        V[intx_data.second], V[neighbors[(n+1)%neighbors.size()]]);
            }
            // TODO: We need to get the incoming angle projected BACK onto the corner of the walked on face,
            // NOT the angle on the Newell face.
            // Start at outgoing edges of the current face corner and add up angles
            for (int n = 0; n < neighbors.size(); n++) {
                k += Utils::vectorAngle(V[intx_data.second], V[neighbors[n]],
                                        V[intx_data.second], V[neighbors[(n+1)%neighbors.size()]]);
            }
            // Get angles beginning from the incoming direction and walk around the fan until we reach halfway.
            // Compute the new direction
            // After attaining new direc, project back onto the tangent plane of the new face (i.e., new face may be nonplanar)

            // TODO: Check if our new direction ends up hitting a boundary. If so, terminate.
            if (-1) {
                return false;
            }

            // Otherwise, use this to recurse
            
            return true;
        } else {    // Landed on edge
            // Query the adjacent faces and grab the one that isn't the current
            std::vector<int> faceNeighbors;
            heMesh.edgeFaces(intx_data.second, faceNeighbors);
            std::pair<int, int> vertNeighbors = heMesh.he_index2directed_edge(associatedHalfEdge);
            // TODO: Compress this if else statement
            if (faceNeighbors[0] == el_idx) {
                new_el_idx = faceNeighbors[1];
            } else {
                new_el_idx = faceNeighbors[0];
            }
            // If the adjacent face is a boundary, terminate
            if (el_idx == -1) {
                return false;
            }

            // Project edge onto new face's tangent plane
            Eigen::Vector3d proj;
            Utils::projectVectorOntoTangentPlane(fNormals[new_el_idx], V[vertNeighbors.first] - V[vertNeighbors.second], proj);
            // Rotate about the start point (vertNeighbors.second) by an angle of theta using the normal.
            // The normalized result is the new direc


            return true;
        }
    }

}   // namespace Mesh