#include "mesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <limits>
#include <utility>

// Utility functions for mesh (projection, insertion, etc.)

namespace Mesh {

// Project a vertex onto the mesh. If multiple, just picks the one with smaller index.
// Also returns the element type that was landed on.
// For non-planar faces, I am just going to fit a Newell plane using the barycenter and vector area + a barycentric height interpolation
int mesh::computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, int& elIdx, bool snap) {
    double tol = 1e-6 * bboxDiag;
    double min_dist = std::numeric_limits<double>::infinity();

    int closest_f = -1;
    elIdx = -1;
    if (active_f == 0) {
        return -1;
    }
    // First find closest face by iterating over faces
    // NOTE: For triangles, this can be done much more simply using
    // barycentric coordinates w/ a linear solve. For arbitrary non-planar polygons,
    // this isn't possible, since polygons may not be convex
    for (int f = 0; f < F.size(); f++) {
        int local_elType = 2;
        int local_elIdx = -1;
        Eigen::Vector3d v_proj;
        // Skip inactive faces
        if (!F[f].active) {
            continue;
        }
        std::vector<int> fVerts = F[f].verts;
        int fSize = fVerts.size();
        Eigen::Vector3d fNormal = F[f].n;
        std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(f);
        // If the face is a triangle, just compute triangle closest point
        if (fVerts.size() == 3) {
            v_proj = Utils::triangleClosestPoint(fVertsPos, v);
        } else {    // Otherwise, we need to use the Newell plane approach
            // Get barycenter 
            Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertsPos);

            // Build local 2D basis
            Eigen::Vector3d t1;
            Eigen::Vector3d t2;
            Utils::buildPlaneBasis(fNormal, t1, t2);
            // Project p onto Newell plane and get its 2D coordinate
            Eigen::Vector3d v_proj3d = Utils::projectPointOntoPlane(fNormal, barycenter, v);
            Eigen::Vector2d v_proj2d = Utils::convertTo2D(v_proj3d, barycenter, t1, t2);

            // Project face vertices onto Newell plane using basis vectors
            std::vector<Eigen::Vector2d> fVert2D(fSize);
            for (int fv = 0; fv < fSize; fv++) {
                Eigen::Vector3d fv_proj3D = Utils::projectPointOntoPlane(fNormal, barycenter, fVertsPos[fv]);
                fVert2D[fv] = Utils::convertTo2D(fv_proj3D, barycenter, t1, t2);
            }
            
            Eigen::Vector2d v_cp;
            // Check if 2D point is in Newell polygon
            // If so, take that one
            if (Utils::pointInPolygon2D(v_proj2d, fVert2D)) {
                v_cp = v_proj2d;
            } else {
                // Project onto all edges to find closest point in 2D
                double min_eDist = std::numeric_limits<double>::infinity();
                for (int i = 0; i < fVert2D.size(); i++) {
                    int j = (i + 1) % fVert2D.size();
                    Eigen::Vector2d cp = Utils::closestPointOnSegment2D(v_proj2d, fVert2D[i], fVert2D[j]);
                    double eDist = (cp - v_proj2d).squaredNorm();
                    if (eDist < min_eDist) {
                        min_eDist = eDist;
                        v_cp = cp;
                    }
                }
            }
            v_proj = Utils::revertTo3D(v_cp, barycenter, t1, t2);
        }
        // True distance in 3D from query point
        double dist = (v - v_proj).norm();
        if (dist < min_dist) {
            min_dist = dist;
            proj = v_proj;
            elIdx = f;
        }
    }

    // Definitive closest face's data
    std::vector<int> fVerts = F[elIdx].verts;
    int fSize = fVerts.size();
    std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(elIdx);
    Eigen::Vector3d fN = F[elIdx].n;

    // Snap to nearby vertex or edge if we are too close
    // NOTE: using Euclidean, not geodesic distance, since this is too hard for non-planar faces
    if (snap) {
        // For the face that was landed on, check if we are close to a vertex on that face
        for (int fv = 0; fv < fSize; fv++) {
            // Once found, we can just return immediately
            if ((proj - fVertsPos[fv]).norm() <= tol) {
                proj = fVertsPos[fv];
                elIdx = fVerts[fv];
                return 0;
            }
        }

        // If not, check if we are close to an edge in 3D
        for (int fv = 0; fv < fSize; fv++) {
            int fv1 = (fv + 1) % fSize;
            Eigen::Vector3d v_projE = Utils::closestPointOnSegment3D(proj, fVertsPos[fv], fVertsPos[fv1]);
            // Once found, we can just return immediately
            if ((proj - v_projE).norm() <= tol) {
                proj = v_projE;
                int he_idx = vertPairToHE.at({fVerts[fv], fVerts[fv1]});
                elIdx = HE[he_idx].edge;
                return 1;
            }
        }
    }

    // If not snapping, then we must be on a face
    // Check if we are on a non-planar face. If so, pin-point the location using MVC
    Eigen::VectorXd fHeight = computeFaceHeight(elIdx);
    bool planar = true;
    for (int v = 0; v < fSize; v++) {
        if (std::abs(fHeight(v)) >= 1e-6) {
            planar = false;
        }
    }
    if (planar || fHeight.size() == 3) {   // Planar face, no MVC interpolation to be done
        return 2;
    }

    // Otherwise, we need to compute mean value coordinates to get projection location
    Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertsPos);
    // Build local 2D basis
    Eigen::Vector3d t1 = (fVertsPos[1]- fVertsPos[0]).normalized();
    Eigen::Vector3d t2;
    Utils::buildPlaneBasis(fN, t1, t2);
    Eigen::Vector2d v_proj2d = Utils::convertTo2D(proj, barycenter, t1, t2);

    std::vector<Eigen::Vector2d> fVerts2D(fSize);
    for (int fv = 0; fv < fSize; fv++) {
        Eigen::Vector3d fv_proj3D = Utils::projectPointOntoPlane(fN, barycenter, fVertsPos[fv]);
        fVerts2D[fv] = Utils::convertTo2D(fv_proj3D, barycenter, t1, t2);
    }
    // Second, apply mean value coordinates to get height function weights
    Eigen::VectorXd MVCWeights(fSize);
    Utils::meanValueCoordinates(v_proj2d, fVerts2D, MVCWeights);

    // Now recover the height using MVC weights
    double h = MVCWeights.dot(fHeight);
    // Add height to current Newell projection
    proj += h * fN;

    return 2;
}


// Inserts a vertex at a location into a data structure
// Returns the index of the new vertex
// If no face information is provided, just -1
// NOTE: The new vertex has no normal information or associated halfedge.
int mesh::insertVertex(Eigen::Vector3d pos, int f, int dCN_idx) {
    int v = V.size();
    V.emplace_back();
    V[v].pos = pos;
    active_v++;
    // Augment the associated face
    if (f != -1) {
        F[f].verts.push_back(v);
    }
    // Add in the dCN info
    if (dCN_idx != -1) {
        V[v].dCN_idx = dCN_idx;
    }
    return v;
}

// Topologically splits an existing edge by adding a new vertex.
// NOTE: Added vertex does NOT need to lie on the edge
// Returns index of the new vertex
int mesh::splitEdge(int e, Eigen::Vector3d split_pos) {
    // Make and modify a copy of the current halfedges
    int he0_idx = E[e].he;
    int he1_idx = HE[he0_idx].twin;
    HalfEdge he0 = HE[he0_idx];
    HalfEdge he1 = HE[he1_idx];
    int u = HE[he1_idx].dest;
    int v = HE[he0_idx].dest;
    // Insert the new vertex into the mesh
    int new_v = insertVertex(split_pos);
    // Insert the new halfedges
    int he0_new = HE.size();
    int he1_new = he0_new + 1;
    HE.emplace_back();
    HE.emplace_back();
    HE[he0_new] = he0;
    HE[he1_new] = he1;
    // Insert the new edge
    int new_e = E.size();
    E.emplace_back();
    E[new_e].he = he0_new;
    active_e++;

    // Re-wire the halfedges
    HE[he0_new].twin = he1_new;
    HE[he1_new].twin = he0_new;
    HE[he0_new].edge = new_e;
    HE[he1_new].edge = new_e;
    HE[he0_idx].dest = new_v;
    HE[he1_new].dest = new_v;

    // Check for how the next and previous halfedges cycle (i.e., how he0_new relates to he1_new)
    bool isEnd = false;
    if (HE[he0_idx].next == HE[he0_idx].twin) {  // Current edge cuts into mesh
        HE[he0_new].next = he1_new;
        HE[he1_new].prev = he0_new;
        isEnd = true;
    }
    if (HE[he1_idx].next == HE[he0_idx].twin) { // Two-sided to allow for splitting standalone edge
        HE[he1_new].next = he0_new;
        HE[he0_new].prev = he1_new;
        isEnd = true;
    }
    if (!isEnd) {    // Standard mesh edge
        HE[HE[he0_idx].next].prev = he0_new;
        HE[HE[he1_idx].next].prev = he1_new;
    }
    // Necessary rewiring between old and new
    HE[he0_idx].next = he0_new;
    HE[he0_new].prev = he0_idx;

    HE[he1_new].next = he1_idx;
    HE[he1_idx].prev = he1_new;

    V[new_v].he = he0_new;
    if (HE[he1_idx].face == -1) {   // Preserve outgoing boundary HE
        V[new_v].he = he1_idx;
    }
    // Modify vertPairToHE
    // Erase old, and insert replacements
    vertPairToHE.erase({u, v});
    vertPairToHE.erase({v, u});
    vertPairToHE[{u, new_v}] = he0_idx;
    vertPairToHE[{new_v, u}] = he1_idx;
    // New keys
    vertPairToHE[{new_v, v}] = he0_new;
    vertPairToHE[{v, new_v}] = he1_new;

    // Update adjacent non-boundary faces' loops and recompute normals and areas
    if (HE[he0_idx].face != -1) {
        F[HE[he0_idx].face].verts.push_back(new_v);
        Eigen::Vector3d n;
        F[HE[he0_idx].face].fArea = computeFVectorArea(HE[he0_idx].face, n);
        F[HE[he0_idx].face].n = n;
    } if (HE[he1_idx].face != -1) {
        if (HE[he0_idx].face != HE[he1_idx].face) { // Do not insert duplicate vertices
            F[HE[he1_idx].face].verts.push_back(new_v);
        }
        Eigen::Vector3d n;
        F[HE[he1_idx].face].fArea = computeFVectorArea(HE[he1_idx].face, n);
        F[HE[he1_idx].face].n = n;
    }

    return new_v;
}

// Inserts a new edge connecting two vertices on a specified face
// Returns the index of the new edge
int insertEdge(int f, int v0, int v1, int dCN_idx, int dCN_idx0, int dCN_idx1, bool positive0) {
    // Check if both vertices actually share the specified face and do not have an existing edge (or if at least one of the vertices has no halfedge)

    // If so, create new halfedges he0 and he1 connecting the two

    // Create a new edge for the halfedges

    // Copy the old next0, prev0; old next1, prev1 and reconnect using the new he0 and he1
    // BE CAREFUL: if next0.twin == prev1 of next1.twin == prev0, then that means we are adding an edge off an endpoint --> CANNOT close the face
    // Another case here?

    // Starting from he0, trace the new halfedge path until either he1 is reached, or he0 is reached.
    // If he0 is reached first, then trace he1 as well. Split into two faces.
    // If he1 is reached first, then we have some sort of endpoint structure or an enclosed shape --> no split
    //      Finish the loop and keep track of the entire loop.
    //      ???
    //      I think checking for enclosed new faces might be too expensive. Warn user not to produce this bad behavior

    // Face splitting
    // 1. Create a new face for h1. For every single HE in this loop, change its face index to the new face index
    //      For the new face index, add in all of the new vertices that are dest points of the HE loop.
    // 2. For the old face, simply remove all of the vertices that are dest points of the new face.

    // If NOT splitting, make sure to add any new vertices to the face's vertex list.

    // return the new edge index
    return 1;
}


}   // namespace Mesh