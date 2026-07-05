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

// Trace a "straightest" geodesic (ish) from the start vert to the end
int mesh::traceGeodesic(const Vert& start, 
                  const Vert& end, 
                  Eigen::Vector3d direc, 
                  int walk_ElType,
                  int walk_ElIdx,
                  std::vector<Vert>& tracedVerts,
                  bool recompute,
                  bool fast) {
    // First, do some simple tests for termination
    // It's good to have these to catch tiny directional drift
    if (walk_ElType == end.mesh_elType && walk_ElIdx == end.mesh_elIdx) {   // The next mesh element is exactly the goal
        return true;
    } else if (start.mesh_elType == 0 && end.mesh_elType == 0) {    // Both are vertices
        if (vertPairToHE.find(std::make_pair(start.mesh_elIdx, end.mesh_elIdx)) != vertPairToHE.end()) {
            return true;
        }
    } else if (start.mesh_elType == 0 && end.mesh_elType == 1) {    // Start is vertex, end is edge
        // Check if either end of the edge is the vertex
        if (HE[E[end.mesh_elIdx].he].dest == start.mesh_elIdx || HE[HE[E[end.mesh_elIdx].he].twin].dest == start.mesh_elIdx) {
            return true;
        }
    } else if (start.mesh_elType == 1 && end.mesh_elType == 0) {    // Start is edge, end is vertex
        // Check if either end of the edge is the vertex
        if (HE[E[start.mesh_elIdx].he].dest == end.mesh_elIdx || HE[HE[E[start.mesh_elIdx].he].twin].dest == end.mesh_elIdx) {
            return true;
        }
    } else if (start.mesh_elType == 1 && end.mesh_elType == 1) {    // Both are on edges
        // Check if they share an edge
        if (start.mesh_elIdx == end.mesh_elIdx) {
            return true;
        }
    }

    // Otherwise need to do a face-wise check
    std::vector<int> startAdjF = adjFaces(start.mesh_elType, start.mesh_elIdx);
    std::vector<int> endAdjF = adjFaces(end.mesh_elType, end.mesh_elIdx);

    std::vector<int> sharedAdjF;
    for (int i = 0; i < startAdjF.size(); i++) {
        if (startAdjF[i] == -1) {
            continue;
        }
        for (int j = 0; j < endAdjF.size(); j++) {
            if (endAdjF[j] == -1) {
                continue;
            }
            if (startAdjF[i] == endAdjF[j]) {
                sharedAdjF.push_back(startAdjF[i]);
            }
        }
    }
    // Check if there are any shared faces (termination condition)
    if (sharedAdjF.size() >= 1) {   // Share at least one face
        if (fast) { // Fast version is a simple face check
            return true;
        } else {    // Slow version does a visibility check
            for (int f_idx : sharedAdjF) {
                if (testVisibility(f_idx, start.pos, end.pos)) {
                    return true;
                }
            }
        }
    }

    // If we failed all those tests, then we actually need to walk! :(
    // First, recompute direction if needed
    if (recompute) {
        direc = (end.pos - start.pos).normalized();
    }

    // Case 1: Check if walk direction is on an edge, simply grab the other end vertex of the edge
    if (walk_ElType == 1) {
        int he = E[walk_ElIdx].he;
        Eigen::Vector3d heVec = (V[HE[he].dest].pos - V[HE[HE[he].twin].dest].pos).normalized();
        Eigen::Vector3d projDirec = direc.dot(heVec) * heVec;
        // Orient ourselves correctly
        int next;
        if (projDirec.dot(heVec) >= 0) {
            next = HE[HE[he].twin].dest;
        } else {
            next = HE[he].dest;
        }
        Vert nextVert = createVertex(V[next].pos, V[next].n, 2, -1, 0, next);
        // Compute the next direction at the intersection

        // TODO: Specially handling for boundaries
        Eigen::Vector3d nextDirec;
        int next_ElType;
        int next_ElIdx;
        // Recurse
        tracedVerts.push_back(nextVert);
        traceGeodesic(nextVert, end, nextDirec, next_ElType, next_ElIdx, tracedVerts, true, fast);
    }

    // If we reached this point, we are definitely walking on a face
    // Project walk direction onto specified direction
    Eigen::Vector3d hit;
    int hit_ElIdx;
    int hit_ElType = rayCastOnFace(walk_ElIdx, start.pos, direc, hit, hit_ElIdx);
    if (hit_ElType == -1) {
        return -1;
    }
    // Produce a "corrected" direction to account for projection drift
    Eigen::Vector3d corrDirec = (hit - start.pos).normalized();
    // Create a new vertex at intersection and append to list
    Vert nextVert = createVertex(hit, getNormal(hit_ElType, hit_ElIdx), 2, -1, hit_ElType, hit_ElIdx);

    int next_ElType, next_ElIdx;
    Eigen::Vector3d nextDirec;
    // TODO: Figure out the next direction based on the intersection. Specially handling for boundaries!

    tracedVerts.push_back(nextVert);
    traceGeodesic(nextVert, end, nextDirec, next_ElType, next_ElIdx, tracedVerts, true, fast);
}

// Test whether the start and end are visible from each other on a particular face
// Returns true if so.
bool mesh::testVisibility(int f, Eigen::Vector3d start, Eigen::Vector3d end, double eps) {
    std::vector<Eigen::Vector3d> fVerts = faceAdjVerts(f);
    if (fVerts.size() == 3) {   // Triangles must be convex
        return true;
    }
    Eigen::Vector3d t1, t2;
    Utils::buildPlaneBasis(F[f].n, t1, t2);
    // For simplicity, assume tangent plane is centered on the start
    vector2dList projFVerts(fVerts.size());
    Eigen::Vector2d start2D = Eigen::Vector2d::Zero();
    Eigen::Vector2d end2D = Utils::convertTo2D(end, start, t1, t2);
    for (int v = 0; v < fVerts.size(); v++) {
        Eigen::Vector3d proj3D = Utils::projectPointOntoPlane(F[f].n, start, fVerts[v]);
        projFVerts[v] = Utils::convertTo2D(proj3D, start, t1, t2);
    }
    Eigen::Vector2d direc = end2D - start2D;
    double t_end = direc.norm();
    direc.normalize();

    std::vector<double> intersections;
    // Raycast to find nearest segment, tracking the t-vals of intersections
    for (int v = 0; v < projFVerts.size(); v++) {
        int v_p1 = (v+projFVerts.size()+1) % projFVerts.size();
        double t, u;
        // Raycast to find nearest segment
        if (Utils::raycastToSegment2D(start2D, direc, projFVerts[v], projFVerts[v_p1], t, u, false)) {
            // Filter out any that are behind
            if (t >= eps) {
                intersections.push_back(t);
            }
        }
    }

    // If intersects none then we are outside the face and/or pointing the wrong way... not good, but check anyways
    if (intersections.size() == 0) {
        return false;
    }
    // Otherwise, check if we hit anything first
    for (double t : intersections) {
        // If we are w/in an epsilon, then ignore
        if (t <= t_end - eps) {
            return false;
        }
    }
    return true;
}

// Find the next intersection point while walking on a particular face
int mesh::rayCastOnFace(int f, Eigen::Vector3d start, Eigen::Vector3d direc, Eigen::Vector3d& hit, int& hit_ElIdx, double eps) {
    double tol = 1e-6 * bboxDiag;
    Eigen::Vector3d projDirec;
    Utils::projectVectorOntoTangentPlane(F[f].n, direc, projDirec);
    std::vector<int> fVertIdxs = faceAdjVertIdxs(f);
    std::vector<Eigen::Vector3d> fVerts = faceAdjVerts(fVertIdxs);

    Eigen::Vector3d t1, t2;
    Utils::buildPlaneBasis(F[f].n, t1, t2);
    // For simplicity, assume tangent plane is centered on the start
    vector2dList projFVerts(fVerts.size());
    Eigen::Vector2d start2D = Eigen::Vector2d::Zero();
    Eigen::Vector2d direc2D = Utils::convertTo2D(projDirec, start, t1, t2);
    for (int v = 0; v < fVerts.size(); v++) {
        Eigen::Vector3d proj3D = Utils::projectPointOntoPlane(F[f].n, start, fVerts[v]);
        projFVerts[v] = Utils::convertTo2D(proj3D, start, t1, t2);
    }
    direc2D.normalize();

    // intersections are a pair of local segment index and the u parameter along that segment
    std::vector<std::pair<int, double>> intersections;
    std::vector<double> intersections_t;
    // Raycast to find nearest segment, tracking the t-vals of intersections
    for (int v = 0; v < projFVerts.size(); v++) {
        int v_p1 = (v+projFVerts.size()+1) % projFVerts.size();
        double t, u;
        // Raycast to find nearest segment
        if (Utils::raycastToSegment2D(start2D, direc2D, projFVerts[v], projFVerts[v_p1], t, u, true)) {
            // Filter out any that are behind
            if (t >= eps) {
                intersections.push_back(std::make_pair(v, u));
                intersections_t.push_back(t);
            }
        }
    }

    // If intersects none then we are outside the face and/or pointing the wrong way... not good, but check anyways
    if (intersections.size() == 0) {
        return -1;
    }
    // Otherwise, check what we hit first
    int nearest_idx = 0;
    double nearest = intersections_t[0];
    for (int i = 0; i < intersections_t.size(); i++) {
        if (intersections_t[i] <= nearest) {
            nearest = intersections_t[i];
            nearest_idx = i;
        }
    }
    // Compute the intersection
    int v = intersections[nearest_idx].first;
    double u = intersections[nearest_idx].second;
    hit = fVerts[v] + u * (fVerts[(v+1)%fVerts.size()] - fVerts[v]);
    // Apply snapping as necessary
    int next_ElType = 1;
    int next_ElIdx = HE[vertPairToHE[std::make_pair(fVertIdxs[v], fVertIdxs[(v+1)%fVertIdxs.size()])]].edge;
    if ((hit - fVerts[v]).norm() <= tol) {
        hit = fVerts[v];
        next_ElType = 0;
        next_ElIdx = fVertIdxs[v];
    }
    return next_ElType;
}

}   // namespace Mesh