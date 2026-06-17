#include "mesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <limits>
#include <utility>
#include <algorithm>

// Utility functions for mesh (projection, insertion, etc.)

namespace Mesh {

// Project a vertex onto the mesh. If multiple, just picks the one with smaller index.
// Also returns the element type that was landed on.
// For non-planar faces, I am just going to fit a Newell plane using the barycenter and vector area + a barycentric height interpolation
int mesh::computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, int& elIdx, bool snap, bool fast) const {
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
        std::vector<int> fVerts = faceAdjVertIdxs(f);
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
            vector2dList fVert2D(fSize);
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
    std::vector<int> fVerts = faceAdjVertIdxs(elIdx);
    int fSize = fVerts.size();
    std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(fVerts);
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
    // If using the fast version, just take the current Newell plane nearest
    if (fast) {
        return 2;
    }
    // Else do the slow way: lift proj using height
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
    Eigen::Vector3d t1 = (fVertsPos[1] - fVertsPos[0]).normalized();
    Eigen::Vector3d t2;
    Utils::buildPlaneBasis(fN, t1, t2);
    Eigen::Vector2d v_proj2d = Utils::convertTo2D(proj, barycenter, t1, t2);

    vector2dList fVerts2D(fSize);
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

// Inserts a vertex at a location and face
// Returns the index of the new vertex
// NOTE: The new vertex has no normal information or associated halfedge.
int mesh::insertVertex(Eigen::Vector3d pos, int f, int dCN_idx) {
    // Augment the associated face
    // Do not accept any invalid or inactive faces
    if (f < 0 || f >= F.size() || !F[f].active) {
        return -1;
    }
    int v = V.size();
    V.emplace_back();
    V[v].pos = pos;
    active_v++;
    F[f].verts.push_back(v);
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
    int f0_idx = HE[he0_idx].face;
    int f1_idx = HE[he1_idx].face;
    HalfEdge he0 = HE[he0_idx];
    HalfEdge he1 = HE[he1_idx];
    int u = HE[he1_idx].dest;
    int v = HE[he0_idx].dest;

    // Insert the new vertex into the mesh
    int new_v;
    if ((f0_idx == f1_idx) && (f0_idx >= 0)) {  // No duplicate insertions
        new_v = insertVertex(split_pos, f0_idx);
    } else if (f0_idx >= 0) {
        new_v = insertVertex(split_pos, f0_idx);
        if (f1_idx >= 0) {
            F[f1_idx].verts.push_back(new_v);
        }
    } else if (f1_idx >= 0) {
        new_v = insertVertex(split_pos, f1_idx);
        if (f0_idx >= 0) {
            F[f0_idx].verts.push_back(new_v);
        }
    } else {
        return -1;  // i.e., not a valid edge to split
    }
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
    // NOTE: Since it's copied from he0 and he1, we don't need to rewire halfedge faces
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

    // Recompute adjacent face normals and areas
    if (f0_idx != -1) {
        Eigen::Vector3d n;
        F[f0_idx].fArea = computeFVectorArea(HE[he0_idx].face, n);
        F[f0_idx].n = n;
    } if (f1_idx != -1) {
        Eigen::Vector3d n;
        F[HE[he1_idx].face].fArea = computeFVectorArea(HE[he1_idx].face, n);
        F[HE[he1_idx].face].n = n;
    }

    return new_v;
}

int mesh::insertEdge(int f, int v0, int v1, int dCN_idx0, int dCN_idx1, bool positive0) {
    // Validate inputs
    if (f < 0 || f >= F.size() || !F[f].active) {   // valid face
        return -1;
    }
    if (v0 < 0 || v0 >= V.size() || !V[v0].active) {    // valid v0
        return -1;
    }
    if (v1 < 0 || v1 >= V.size() || !V[v1].active) {    // valid v1
        return -1;
    }
    if (v0 == v1) { // cannot draw edge from a point to itself
        return -1;
    }

    // Check if edge already exists
    if (vertPairToHE.find({v0, v1}) != vertPairToHE.end() || vertPairToHE.find({v1, v0}) != vertPairToHE.end()) {
        return -1;
    }

    // Check that both verts belong to this face
    bool v0_inF = false;
    bool v1_inF = false;
    for (int fv : F[f].verts) {
        if (fv == v0) {
            v0_inF = true;
        }
        if (fv == v1) {
            v1_inF = true;
        }
    }
    if (!v0_inF || !v1_inF) {   // If they don't belong to the face, then return
        return -1;
    }

    // Note if any vertices are isolated
    const bool v0_isolated = (V[v0].he == -1);
    const bool v1_isolated = (V[v1].he == -1);  
    int he0_prev = -1;
    int he1_prev = -1;

    // ------------------------------------------------------------
    // Slot selection
    // ------------------------------------------------------------
    if (!v0_isolated) {
        if (!chooseEdgeInsertHE(f, v0, v1, he0_prev)) {
            return -1;
        }
    }
    if (!v1_isolated) {
        if (!chooseEdgeInsertHE(f, v1, v0, he1_prev)) {
            return -1;
        }
    }

    // Based on these, apply edge insertion
    return insertEdgeBetweenHEs(
        f,
        v0,
        v1,
        he0_prev,
        v0_isolated,
        he1_prev,
        v1_isolated,
        dCN_idx0,
        dCN_idx1,
        positive0
    );
}

// Inserts a new edge connecting two vertices on a specified face
// Returns the index of the new edge
// NOTE: This should only mainly be used as a helper for insertEdge()
// NOTE: Paper says not to update normals here. Still included, but maybe try taking out later?
int mesh::insertEdgeBetweenHEs(int f, int v0, int v1, int he0_prev, int v0_isolated, int he1_prev, int v1_isolated, int dCN_idx0, int dCN_idx1, bool positive0) {
    // Create new halfedges he0 and he1 connecting the two vertices
    int he0 = HE.size();
    int he1 = he0+1;
    HE.emplace_back();
    HE.emplace_back();
    // Convention: u->v is the pos/0 side, and v->u is the neg/1 side
    HE[he0].dest = v1;
    HE[he1].dest = v0;
    HE[he0].twin = he1;
    HE[he1].twin = he0;
    // Replace vertex halfedges if not already pointing to a boundary
    if (v0_isolated || !vertIsBoundary(v0)) {
        V[v0].he = he0;
    }
    if (v1_isolated || !vertIsBoundary(v1)) {
        V[v1].he = he1;
    }
    // For safety, set initial face to f
    HE[he0].face = f;
    HE[he1].face = f;
    // Separate the index assignment of dCN for both halfedges
    if (dCN_idx0 != -1) {
        HE[he0].dCN_idx = dCN_idx0;
        HE[he0].dCN_sign = positive0;
    }
    if (dCN_idx1 != -1) {
        HE[he1].dCN_idx = dCN_idx1;
        HE[he1].dCN_sign = !positive0;
    }
    vertPairToHE[{v0, v1}] = he0;
    vertPairToHE[{v1, v0}] = he1;
    // Create a new edge for the halfedges
    int e = E.size();
    E.emplace_back();
    E[e].he = he0;
    HE[he0].edge = e;
    HE[he1].edge = e;
    active_e++;

    // Case 1: both vertices are isolated
    if (v0_isolated && v1_isolated) {
        HE[he0].next = he1;
        HE[he0].prev = he1;
        HE[he1].next = he0;
        HE[he1].prev = he0;
        return e;
    }

    // NOTE: The following if statements can likely be streamlined
    // Case 2: v0 is isolated but v1 is not
    else if (v0_isolated && !v1_isolated) {
        HE[he0].next = HE[he1_prev].next;
        HE[HE[he1_prev].next].prev = he0;

        HE[he0].prev = he1;
        HE[he1].next = he0;

        HE[he1].prev = he1_prev;
        HE[he1_prev].next = he1;
        return e;
    } else if (!v0_isolated && v1_isolated) {
        HE[he1].next = HE[he0_prev].next;
        HE[HE[he0_prev].next].prev = he1;

        HE[he0].next = he1;
        HE[he1].prev = he0;

        HE[he0].prev = he0_prev;
        HE[he0_prev].next = he0;
        return e;
    } else {
        // Otherwise, both are not isolated
        HE[he0].next = HE[he1_prev].next;
        HE[he0].prev = he0_prev;
        HE[he1].next = HE[he0_prev].next;
        HE[he1].prev = he1_prev;

        HE[he0_prev].next = he0;
        HE[HE[he1].next].prev = he1;
        HE[he1_prev].next = he1;
        HE[HE[he0].next].prev = he0;

        // Check if we need to split the face
        // If we trace a loop from he0 and can recover he1, then we must form a closed face
        // Otherwise, the two must be disconnected
        bool split_face = true;
        std::vector<int> f0_HEloop = halfedgeLoop(he0);
        for (int i = 0; i < f0_HEloop.size(); i++) {
            if (f0_HEloop[i] == he1) {
                split_face = false;
                break;
            }
        }
        if (split_face) {
            // Face Split strategy: Let he0 and its loop keep the original face
            // Create a new face for he1's loop. Then remove the vertices in the new face from the vert list of the original.
            int new_f = F.size();
            F.emplace_back();
            F[new_f].he = he1;  // Always give the new face to the cycle containing he1
            HE[he1].face = new_f;

            std::vector<int> f1_HEloop = halfedgeLoop(he1);
            // Get new face's vertices
            std::vector<int> f1_Vloop(f1_HEloop.size());
            for (int i = 0; i < f1_HEloop.size(); i++) {
                f1_Vloop[i] = HE[f1_HEloop[i]].dest;
            }
            // Remove unique vertices in new face from the first face
            for (int f1_he = 0; f1_he < f1_HEloop.size(); f1_he++) {
                HE[f1_HEloop[f1_he]].face = new_f;
                int he_vert = HE[f1_HEloop[f1_he]].dest;
                // Check if the vertex is shared by f0 (i.e., twin's face)
                if (HE[HE[f1_HEloop[f1_he]].twin].face == f) {
                    continue;   // Skip, do not delete
                }
                int rm_idx = -1;
                // Find removal index
                for (int f0_v = 0; f0_v < F[f].verts.size(); f0_v++) {
                    if (he_vert == F[f].verts[f0_v]) { // Check if the vertex list in f contains f1_v
                        rm_idx = f0_v;
                    }
                }
                if (rm_idx != -1) { // Remove
                    F[f].verts.erase(F[f].verts.begin() + rm_idx);
                }
            }
            // Assign the new face its vertex loop
            F[new_f].verts = f1_Vloop;
            F[new_f].he = he1;
            F[f].he = he0;
            active_f++;
            // Recompute face normal and area for adjacent faces
            Eigen::Vector3d f0_n, f1_n;
            double f0_area = computeFVectorArea(f, f0_n);
            double f1_area = computeFVectorArea(new_f, f1_n);
            F[f].fArea = f0_area;
            F[f].n = f0_n;
            F[new_f].fArea = f1_area;
            F[new_f].n = f1_n;
            return e;
        }

        // If NOT splitting, make sure to add any new vertices to the face's vertex list.
        for (int i = 0; i < f0_HEloop.size(); i++) {
            int he_vert = HE[f0_HEloop[i]].dest;
            bool exists = false;
            // Check if vert already exists in the current face's list.
            for (int f0_v = 0; f0_v < F[f].verts.size(); f0_v++) {
                if (he_vert == F[f].verts[f0_v]) { // Check if the vertex list in f contains the vert
                    exists = true;
                }
            }
            if (!exists) {
                F[f].verts.push_back(he_vert);
            }
        }
    }
    // Recompute face normal and area
    Eigen::Vector3d f_n;
    double f_area = computeFVectorArea(f, f_n);
    F[f].fArea = f_area;
    F[f].n = f_n;

    // return the new edge index
    return e;
}

// Helper for edge insertion. Chooses best insertion point
bool mesh::chooseEdgeInsertHE(int f, int v, int target, int& he_prev_out) {
    he_prev_out = -1;

    // Basic validation
    if (f < 0 || f >= F.size() || !F[f].active) {
        return false;
    }
    if (v < 0 || v >= V.size() || !V[v].active) {
        return false;
    }
    if (target < 0 || target >= V.size() || !V[target].active) {
        return false;
    }
    // If a vertex is isolated, immediately return nothing.
    if (V[v].he == -1) {
        return false;
    }

    // Gather outgoing halfedges from v that are incident to face f.
    std::vector<int> adjHEs = vertAdjHEs(v);
    std::vector<int> faceOutgoingHEs;
    for (int he : adjHEs) {
        if (!HE[he].active) {
            continue;
        }
        if (HE[he].face == f) {
            faceOutgoingHEs.push_back(he);
        }
    }
    if (faceOutgoingHEs.empty()) {
        return false;
    }

    // If there's only one incident halfedge to this face, then just grab the previous and return
    if (faceOutgoingHEs.size() == 1) {
        int nextHE = faceOutgoingHEs[0];
        int prevHE = HE[nextHE].prev;
        he_prev_out = prevHE;
        return true;
    }

    // OTHERWISE: We need to do an angle check to figure out
    // where to insert the new halfedge
    const double eps = 1e-12;
    Eigen::Vector3d n = F[f].n;
    // Build 2D basis for Newell plane
    Eigen::Vector3d t1, t2;
    Utils::buildPlaneBasis(n, t1, t2);
    const Eigen::Vector3d& p = V[v].pos;
    // Compute angles on the plane relative to bases and store as a list
    double thetaNew = 0.0;
    if (!Utils::directionAngleInPlane(p, V[target].pos, n, t1, t2, thetaNew)) {
        return false;
    }
    std::vector<std::pair<double, int>> candidateAngles;
    candidateAngles.reserve(faceOutgoingHEs.size());

    for (int he : faceOutgoingHEs) {
        int dst = HE[he].dest;
        double theta = 0.0;
        if (!Utils::directionAngleInPlane(p, V[dst].pos, n, t1, t2, theta)) {
            return false;
        }
        if (Utils::anglesCoincident(theta, thetaNew)) {
            return false;
        }
        candidateAngles.push_back({theta, he});
    }
    // Sort the list by angle
    std::sort(candidateAngles.begin(), candidateAngles.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {return a.first < b.first;});

    // Find first existing outgoing halfedge CCW after the new direction.
    int nextOutgoing = candidateAngles[0].second;
    for (const std::pair<double, int> candidate : candidateAngles) {
        if (candidate.first > thetaNew) {
            nextOutgoing = candidate.second;
            break;
        }
    }
    // Grab the previous halfedge as the one we want
    he_prev_out = HE[nextOutgoing].prev;
    return true;
}

}   // namespace Mesh