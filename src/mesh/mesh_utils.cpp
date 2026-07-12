#include "mesh.hpp"

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
vertProjData mesh::computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, bool snap, bool fast) const {
    double tol = 1e-6 * bboxDiag;
    double min_dist = std::numeric_limits<double>::infinity();
    vertProjData projData({-1, -1});
    int closest_f = -1;
    if (active_f == 0) {
        return projData;
    }
    // First find closest face by iterating over faces
    // NOTE: For triangles, this can be done much more simply using
    // barycentric coordinates w/ a linear solve. For arbitrary non-planar polygons,
    // this isn't possible, since polygons may not be convex
    for (int f = 0; f < F.size(); f++) {
        vertProjData localProjData({2, -1});
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
            projData.elIdx = f;
        }
    }

    // Definitive closest face's data
    std::vector<int> fVerts = faceAdjVertIdxs(projData.elIdx);
    int fSize = fVerts.size();
    std::vector<Eigen::Vector3d> fVertsPos = adjVerts(fVerts);
    Eigen::Vector3d fN = F[projData.elIdx].n;

    // Snap to nearby vertex or edge if we are too close
    // NOTE: using Euclidean, not geodesic distance, since this is too hard for non-planar faces
    if (snap) {
        // For the face that was landed on, check if we are close to a vertex on that face
        for (int fv = 0; fv < fSize; fv++) {
            // Once found, we can just return immediately
            if ((proj - fVertsPos[fv]).norm() <= tol) {
                proj = fVertsPos[fv];
                projData.elIdx = fVerts[fv];
                projData.elType = 0;
                return projData;
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
                projData.elIdx = HE[he_idx].edge;
                projData.elType = 1;
                return projData;
            }
        }
    }

    // If not snapping, then we must be on a face
    // If using the fast version, just take the current Newell plane nearest
    if (fast) {
        return projData;
    }
    // Else do the slow way: lift proj using height
    // Check if we are on a non-planar face. If so, pin-point the location using MVC
    Eigen::VectorXd fHeight = computeFaceHeight(projData.elIdx);
    bool planar = true;
    for (int v = 0; v < fSize; v++) {
        if (std::abs(fHeight(v)) >= 1e-6) {
            planar = false;
        }
    }
    if (planar || fHeight.size() == 3) {   // Planar face, no MVC interpolation to be done
        return projData;
    }

    // Otherwise, we need to compute mean value coordinates to get projection location
    Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertsPos);
    // Build local 2D basis
    Eigen::Vector3d t1 = (fVertsPos[1] - fVertsPos[0]).normalized();
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

    return projData;
}

// Trace a "straightest" geodesic (ish) from the start vert to the end
int mesh::traceGeodesic(const Vert& start, 
                  const Vert& end, 
                  Eigen::Vector3d prevDirec, 
                  int prev_ElType,
                  int prev_ElIdx,
                  std::vector<Vert>& tracedVerts,
                  bool recompute,
                  bool fast) {
    double eps = 1e-6;
    // First, do some simple tests for termination
    // It's good to have these to catch tiny directional drift
    if (prev_ElType == end.projData.elType && prev_ElIdx == end.projData.elIdx) {   // The next mesh element is exactly the goal
        return true;
    } else if (start.projData.elType == 0 && end.projData.elType == 0) {    // Both are vertices
        if (vertPairToHE.find(std::make_pair(start.projData.elIdx, end.projData.elIdx)) != vertPairToHE.end()) {
            return true;
        }
    } else if (start.projData.elType == 0 && end.projData.elType == 1) {    // Start is vertex, end is edge
        // Check if either end of the edge is the vertex
        if (HE[E[end.projData.elIdx].he].dest == start.projData.elIdx || HE[HE[E[end.projData.elIdx].he].twin].dest == start.projData.elIdx) {
            return true;
        }
    } else if (start.projData.elType == 1 && end.projData.elType == 0) {    // Start is edge, end is vertex
        // Check if either end of the edge is the vertex
        if (HE[E[start.projData.elIdx].he].dest == end.projData.elIdx || HE[HE[E[start.projData.elIdx].he].twin].dest == end.projData.elIdx) {
            return true;
        }
    } else if (start.projData.elType == 1 && end.projData.elType == 1) {    // Both are on edges
        // Check if they share an edge
        if (start.projData.elIdx == end.projData.elIdx) {
            return true;
        }
    }

    // Otherwise need to do a face-wise check
    // Find shared faces between start and end, if any
    std::vector<int> startAdjF = adjFaces(start.projData.elType, start.projData.elIdx);
    std::vector<int> endAdjF = adjFaces(end.projData.elType, end.projData.elIdx);
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

    // We have to walk, so we now compute the next walk direction
    int next_ElType, next_ElIdx;
    Eigen::Vector3d nextDirec;
    prevDirec.normalize();
    if (recompute && start.projData.elType == 2) {
        Eigen::Vector3d refDirec = (end.pos - start.pos);
        next_ElType = prev_ElType;
        next_ElIdx = prev_ElIdx;
        // If projection is degenerate, defer to valid previous version
        if (Utils::projectVectorOntoTangentPlane(F[next_ElIdx].n, refDirec, nextDirec) <= eps) {
            nextDirec = prevDirec;
        } else if (nextDirec.dot(prevDirec) < 0.0) {    // Else do a soft check to make sure we're going the right way
            nextDirec *= -1;
        }
    } else {
        if (start.projData.elType == 0) {
            next_ElType = nextEl_Vert(start.projData.elIdx, prev_ElType, prev_ElIdx, prevDirec, nextDirec, next_ElIdx, true);
        } else if (start.projData.elType == 1) {
            next_ElType = nextEl_Edge(start.projData.elIdx, prev_ElIdx, prevDirec, nextDirec, next_ElIdx, true);
        } else {
            next_ElType = prev_ElType;
            next_ElIdx = prev_ElIdx;
            Utils::projectVectorOntoTangentPlane(F[next_ElIdx].n, prevDirec, nextDirec);
        }
    }
    // Screen out any immediately problematic vectors
    if (nextDirec.norm() <= eps) {
        return -1;
    }

    // Case 1: Check if walk direction is on an edge, simply grab the other end vertex of the edge
    if (prev_ElType == 1) {
        // Opt to actually compute the matching direction instead of assuming the start is at a vertex
        // This way we can handle degenerate cases where the initial walk direction is on an edge
        int he = E[prev_ElIdx].he;
        Eigen::Vector3d heVec = (V[HE[he].dest].pos - V[HE[HE[he].twin].dest].pos).normalized();
        Eigen::Vector3d projDirec = (nextDirec.dot(heVec) * heVec).normalized();
        // Orient ourselves correctly
        int next;
        if (projDirec.dot(heVec) >= 0) {
            next = HE[HE[he].twin].dest;
        } else {
            next = HE[he].dest;
        }
        Vert nextVert = createVertex(V[next].pos, V[next].n, 2, -1, 0, next);
        // Recurse
        tracedVerts.push_back(nextVert);
        traceGeodesic(nextVert, end, nextDirec, next_ElType, next_ElIdx, tracedVerts, true, fast);
    }

    // If we reached this point, we are definitely walking on a face
    // Project walk direction onto specified direction
    Eigen::Vector3d hit;
    vertProjData hit_Data;
    int hit_ElIdx;
    int valid = rayCastOnFace(next_ElIdx, start.pos, nextDirec, hit, hit_Data);
    if (valid == -1) {
        return -1;
    }
    // Create a new vertex at intersection and append to list
    Vert nextVert = createVertex(hit, getNormal(hit_Data), 2, -1, hit_Data);

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
    std::vector<Eigen::Vector2d> projFVerts(fVerts.size());
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
        if (Utils::raycastToSegment2D(start2D, direc, projFVerts[v], projFVerts[v_p1], t, u)) {
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
int mesh::rayCastOnFace(int f, Eigen::Vector3d start, Eigen::Vector3d direc, Eigen::Vector3d& hit, vertProjData& hitData, double eps) {
    double tol = 1e-6 * bboxDiag;
    Eigen::Vector3d projDirec;
    Utils::projectVectorOntoTangentPlane(F[f].n, direc, projDirec);
    std::vector<int> fVertIdxs = faceAdjVertIdxs(f);
    std::vector<Eigen::Vector3d> fVerts = adjVerts(fVertIdxs);

    Eigen::Vector3d t1, t2;
    Utils::buildPlaneBasis(F[f].n, t1, t2);
    // For simplicity, assume tangent plane is centered on the start
    std::vector<Eigen::Vector2d> projFVerts(fVerts.size());
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
        if (Utils::raycastToSegment2D(start2D, direc2D, projFVerts[v], projFVerts[v_p1], t, u)) {
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
    for (int i = 1; i < intersections_t.size(); i++) {
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
    hitData.elType = 1;
    hitData.elIdx = HE[vertPairToHE[std::make_pair(fVertIdxs[v], fVertIdxs[(v+1)%fVertIdxs.size()])]].edge;
    if ((hit - fVerts[v]).norm() <= tol) {
        hit = fVerts[v];
        hitData.elType = 0;
        hitData.elIdx = fVertIdxs[v];
    }
    return 1;
}

// Compute the next walk element given that we intersected with an edge
// Returns the next 
int mesh::nextEl_Edge(int e, int f_origin, const Eigen::Vector3d& prev_direc, 
                    Eigen::Vector3d& next_direc, int& next_elIdx, bool bdy_snap) {
    // Get adjacent faces
    std::vector<int> adjF = edgeAdjFaces(e);
    next_elIdx = adjF[0];
    if (f_origin == adjF[0]) {
        next_elIdx = adjF[1];
    }
    // Projection step for safety
    Eigen::Vector3d proj_direc;
    double valid = Utils::projectVectorOntoTangentPlane(F[f_origin].n, prev_direc, proj_direc);
    if (valid <= 0.0) {
        return -1;
    }
    // Handle boundary case first
    // Snap to the better boundary edge direction
    if (next_elIdx == -1) {
        if (bdy_snap) {
            next_elIdx = e;
            // Find the better fit direction 
            int v1 = HE[E[e].he].dest;
            int v0 = HE[HE[E[e].he].twin].dest;
            if (((V[v1].pos - V[v0].pos).normalized()).dot(proj_direc) >= 0) {
                next_direc = (V[v1].pos - V[v0].pos).normalized();
            } else {
                next_direc = (V[v0].pos - V[v1].pos).normalized();
            }
            return 1;
        } else {    // Hit a dead end
            return -1;
        }
    }
    // Otherwise, we just rotate onto the plane of the adajcent face
    // Get the halfedge of the incoming face
    int he0 = E[e].he;
    if (HE[he0].face == f_origin) {
        he0 = HE[he0].twin;
    }
    // Compute the local tangent basis of the origin face
    Eigen::Vector3d tangent = ((V[HE[he0].dest].pos - V[HE[HE[he0].twin].dest].pos).normalized());
    // Normalized for numerical safety
    Eigen::Vector3d biN_origin = tangent.cross(F[f_origin].n).normalized();
    Eigen::Vector3d biN_next = tangent.cross(F[next_elIdx].n).normalized();
    Eigen::Matrix3d rot = Utils::computeRotation(biN_origin, biN_next);
    next_direc = rot * proj_direc;
    return 2;
}

// Compute the next walk element given that we intersected with a vertex
int mesh::nextEl_Vert(int v, int origin_ElType, int origin_ElIdx, 
                    const Eigen::Vector3d& prev_direc, Eigen::Vector3d& next_direc, 
                    int& next_elIdx, bool bdy_snap, double eps) {
    // Gather all of the adjacent halfedges and faces
    std::vector<int> adjHE = vertAdjHEs(v);
    // Figure out which halfedge "matches" the face origin and also find the local umbrella angle sum
    int localIdx = -1;
    double angleSum = 0.0;
    std::vector<double> angleBuckets(adjHE.size());
    std::vector<Eigen::Vector3d> cornerNormals(adjHE.size());
    for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
        int he0 = adjHE[he_idx];
        int he1 = adjHE[(he_idx + 1) % adjHE.size()];
        // Compute angle between the two vectors by first computing an axis
        Eigen::Vector3d vec0 = V[HE[he0].dest].pos - V[v].pos;
        Eigen::Vector3d vec1 = V[HE[he1].dest].pos - V[v].pos;
        Eigen::Vector3d axis = vec0.cross(vec1);
        if (axis.norm() <= eps) {
            axis = F[HE[he0].face].n;
        }
        cornerNormals[he_idx] = axis.normalized();
        // Get the positive signed angle
        angleBuckets[he_idx] = Utils::signedAngle(vec0, vec1, axis, true);
        angleSum += angleBuckets[he_idx];
        // Get the matching face as a local halfedge index in adjHE
        if (origin_ElType == 1 && HE[adjHE[he_idx]].edge == origin_ElIdx) {
            localIdx = he_idx;
            break;
        } else if (origin_ElType == 2 && HE[adjHE[he_idx]].face == origin_ElIdx) {
            localIdx = he_idx;
            break;
        }
    }
    if (localIdx == -1) {   // Did not find a match
        return -1;
    }
    // Compute the projected walk direction onto the corner normal's plane
    Eigen::Vector3d proj_direc;
    double valid = Utils::projectVectorOntoTangentPlane(cornerNormals[localIdx], prev_direc, proj_direc);
    if (valid <= 0.0) {
        return -1;
    }
    double angleToFace = Utils::signedAngle(V[HE[adjHE[localIdx]].dest].pos - V[v].pos, -1 * prev_direc, cornerNormals[localIdx], true);

    // Get the halfway angle
    double targetAngle = angleSum / 2.0;
    int stopHE = -1;
    targetAngle += angleToFace; // Add this on so we can start accumulating angles from the halfedge
    // Iterate over faces to get the target angle
    for (int f = 0; f < adjHE.size(); f++) {
        int f_curr = (localIdx + f) % adjHE.size();
        int f_next = (f_curr + 1) % adjHE.size();
        if (targetAngle <= angleBuckets[f_curr]) {  // Found the face
            stopHE = f_curr;
            break;
        } else {    // Keep going
            targetAngle -= angleBuckets[f_curr];
        }
    }
    // Use the computed face and remnant angle to determine the vector
    Eigen::Vector3d heVec = (V[HE[adjHE[stopHE]].dest].pos - V[v].pos).normalized();
    // To make this well-posed, compute corner normal as the average of adjacent corner normals
    if (cornerNormals[stopHE].norm() <= eps) {
        cornerNormals[stopHE] = (cornerNormals[(stopHE + 1)%adjHE.size()] + cornerNormals[(stopHE + adjHE.size() - 1)%adjHE.size()]) / 2;
    }
    Eigen::Matrix3d rot = Utils::computeRotation(cornerNormals[stopHE], targetAngle);
    next_direc = rot * heVec;

    Eigen::Vector3d heVec1 = (V[HE[adjHE[(stopHE+1)%adjHE.size()]].dest].pos - V[v].pos).normalized();
    next_elIdx = HE[adjHE[stopHE]].face;
    // Process boundary if we hit one
    if (next_elIdx == -1) {
        if (bdy_snap) {
            // Compute best adjacent edge to snap to
            if (next_direc.dot(heVec) <= next_direc.dot(heVec1)) {
                next_direc = heVec;
                next_elIdx = HE[adjHE[stopHE]].edge;
            } else {
                next_direc = heVec1;
                next_elIdx = HE[adjHE[(stopHE+1)%adjHE.size()]].edge;
            }
            return 1;
        } else {    // Failure
            return -1;
        }
    }
    // Check for edge snaps
    if (next_direc.dot(heVec) <= eps) {
        next_direc = heVec;
        next_elIdx = adjHE[stopHE];
        return 1;
    } else if (next_direc.dot(heVec1) <= eps) {
        next_direc = heVec1;
        next_elIdx = adjHE[(stopHE+1)%adjHE.size()];
        return 1;
    }
    return 2;
}

}   // namespace Mesh