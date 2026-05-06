#include "mesh.hpp"

#include "dcurvenet/pdcurvenet.hpp"
#include "polyHE/polyHE.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <limits>
#include <utility>


namespace Mesh {

// Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
mesh::mesh(std::vector<Eigen::Vector3d>& V, std::vector<std::vector<int>>& F): V(V), F(F), num_f(F.size()), num_v(V.size()) {
    heMesh.build_from_face_list(V.size(), F);
    computeFNormalsAreas();
    computeENormals();
    computeVNormalsAreas();
    computeHeightFuncsAndConvexity();
    computeMeanE();
    computeBBoxDiag();
}

// Getters
Eigen::Vector3d& mesh::getVNormal(int vidx) {
    return vNormals[vidx];
}
Eigen::Vector3d& mesh::getENormal(int eidx) {
    return eNormals[eidx];
}
Eigen::Vector3d& mesh::getFNormal(int fidx) {
    return fNormals[fidx];
}
// Get pointer to he mesh
polyHE::polyHE_t& mesh::getHEMesh() {
    return heMesh;
}

// Get mean edge length
double mesh::getMeanE() {
    return meanE;
}

 // Get bbox diagonal
double mesh::getBBoxDiag() {
    return bboxDiag;
}

// Internal function to compute height functions and convexity
void mesh::computeHeightFuncsAndConvexity() {
    // Iterate over faces 
    H.clear();
    H.resize(num_f);
    Convex.clear();
    Convex.resize(num_f);
    for (int f = 0; f < num_f; f++) {
        Convex[f] = true;
        // Special handling for triangles
        if (F[f].size() == 3) {
            H[f] = Eigen::VectorXd({0.0, 0.0, 0.0});
            continue;
        }
        std::vector<Eigen::Vector3d> face_v(F[f].size());
        // Fit Newell Plane. First grab all vertices of face
        for (int v = 0; v < F[f].size(); v++) {
            face_v[v] = V[F[f][v]];
        }
        // 1. compute barycenter and face normal
        Eigen::Vector3d faceCenter = DECUtils::computeBarycenter(face_v);
        Eigen::Vector3d faceN = fNormals[f];
        // 2. Project face vertices onto the Newell plane and grab height
        Eigen::VectorXd faceH = Eigen::VectorXd::Zero(F[f].size());
        std::vector<Eigen::Vector2d> proj_v(F[f].size());
        // Build a basis
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        Utils::buildPlaneBasis(faceN, t1, t2);
        for (int v = 0; v < face_v.size(); v++) {
            Eigen::Vector3d proj3d = Utils::projectPointOntoPlane(faceN, faceCenter, face_v[v]);
            // vector from old point to plane point
            Eigen::Vector3d heightVec = face_v[v] - proj3d;
            // Get 2D version
            proj_v[v] = Utils::convertTo2D(proj3d, faceCenter, t1, t2);
            double height = heightVec.norm();   // How far we are from the plane
            if (height <= 1e-6) {   // If we are on/close to the surface, just snap to the plane
                faceH(v) = 0.0;
            } else if ((heightVec.normalized()).dot(faceN) > 0.0) { // We are above the plane
                faceH(v) = height;
            } else {    // We are below the plane
                faceH(v) = -1.0 * height;
            }
        }
        // Compute signed angles
        for (int v = 0; v < face_v.size(); v++) {
            int v_next = (v+1)%face_v.size();
            int v_prev = (v+face_v.size()-1)%face_v.size();
            double signedAngle = Utils::vectorAngle(proj_v[v], proj_v[v_next], proj_v[v], proj_v[v_prev]);
            // If signed angle is negative AND the signed angle is not close to 180, then probably non-convex
            if ((signedAngle <= 0.0) && (M_PI - signedAngle > 1e-6)) {
                Convex[f] = true;
            }
        }
    }
    return;
}

// Internal function to precompute normals on all mesh structures
void mesh::computeFNormalsAreas() {
    fNormals.clear();
    fNormals.resize(num_f);
    fAreas.clear();
    fAreas.resize(num_f);
    for (int f = 0; f < num_f; f++) {
        // create vector of vertices
        std::vector<Eigen::Vector3d> fList;
        for (int v = 0; v < F[f].size(); v++) {
            fList.push_back(V[F[f][v]]);
        }
        Eigen::Vector3d fN = Eigen::Vector3d::Zero();
        fAreas[f] = DECUtils::vectorArea(fList, fN);
        fNormals[f] = fN;
    }
    return;
}
void mesh::computeENormals(bool weight_fN) {
    eNormals.clear();
    eNormals.resize(heMesh.numEdges());
    for (int e = 0; e < heMesh.numEdges(); ++e) {
        std::vector<int> eFaces;
        heMesh.edgeFaces(e, eFaces);
        int f0 = eFaces[0];
        int f1 = eFaces[1];
        Eigen::Vector3d eN = Eigen::Vector3d::Zero();
        if (weight_fN) {
            if (f0 != -1) {
                eN += fAreas[f0] * fNormals[f0];
            } 
            if (f1 != -1) {
                eN += fAreas[f1] * fNormals[f1];
            }
        } else {
            if (f0 != -1) {
                eN += fNormals[f0];
            }
            if (f1 != -1) {
                eN += fNormals[f1];
            }
        }
        eN.normalize();
        eNormals[e] = eN;
    }
    return;
}
void mesh::computeVNormalsAreas(bool weight_fN) {
    vNormals.clear();
    vNormals.resize(num_v);
    vAreas.clear();
    vAreas.resize(num_v);
    int num_v = V.size();
    for (int v = 0; v < num_v; v++) {
        // create vector of vertices
        std::vector<int> fList;
        heMesh.vertex_face_neighbors(v, fList);
        Eigen::Vector3d vN = Eigen::Vector3d::Zero();
        double vArea = 0.0;
        // Iterate over face list and accumulate areas and normals
        for (int f = 0; f < fList.size(); f++) {
            int fi = fList[f];
            double fArea = fAreas[fi]/(F[fi].size());
            if (weight_fN) {
                vN += fArea * fNormals[fi];
            } else {
                vN += fNormals[fi];
            }
            vArea += fArea;
        }
        vN.normalize();
        vNormals[v] = vN;
        vAreas[v] = vArea;
    }
    return;
}

// Computes mean edge length on the mesh
void mesh::computeMeanE() {
    meanE = 0.0;
    for (int e = 0; e < heMesh.numEdges(); ++e) {
        std::vector<int> eVerts;
        heMesh.edgeVertices(e, eVerts);
        int v0 = eVerts[0];
        int v1 = eVerts[1];
        meanE += (V[v0] - V[v1]).norm();
    }
    meanE /= heMesh.numEdges();
    return;
}

// Computes the diagonal length of the mesh's AABB
void mesh::computeBBoxDiag() {
    Eigen::Vector3d minV = V[0];
    Eigen::Vector3d maxV = V[0];
    // Find most extreme points in mesh
    for (int v = 1; v < V.size(); v++) {
        minV = minV.cwiseMin(V[v]);
        maxV = maxV.cwiseMax(V[v]);
    }
    // get norm of the most extreme points
    bboxDiag = (maxV - minV).norm();
    return;
}

// Given a point in space, computes the index of the nearest edge to that point
double mesh::computeNearestFaceEdge(const int face, const Eigen::Vector3d& p, int& nearestIdx, Eigen::Vector3d& nearestPnt) {
    double min_dist = std::numeric_limits<double>::infinity();
    nearestIdx = -1;
    // Iterate over edges
    std::vector<int> fVList = F[face];
    for (int v = 0; v < fVList.size(); v++) {
        Eigen::Vector3d proj = Utils::closestPointOnSegment3D(p, V[fVList[v]], V[fVList[(v+1)%fVList.size()]]);
        double dist = (proj - p).norm();
        if (dist < min_dist) {
            min_dist = dist;
            nearestIdx = heMesh.edgeIdxFromVerts(fVList[v], fVList[v+1]%fVList.size());
            nearestPnt = proj;
        }
    }
    return min_dist;
}

// Project a vertex onto the mesh. If multiple, just picks the first one.
// Also returns the element type that was landed on
// For simplicity, I am just going to fit a Newell plane using the barycenter and vector area
int mesh::computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, int& elIdx, bool snap) {
    double tol = 1e-6 * bboxDiag;
    double min_dist = std::numeric_limits<double>::infinity();

    int closest_f = -1;
    elIdx = -1;
    // First find closest face by iterating over faces
    // NOTE: For triangles, this can be done much more simply using
    // barycentric coordinates w/ a linear solve. For arbitrary non-planar polygons,
    // this isn't possible, since polygons may not be convex
    for (int f = 0; f < F.size(); f++) {
        int local_elType = 2;
        int local_elIdx = -1;
        
        // Get barycenter 
        std::vector<Eigen::Vector3d> fVertList(F[f].size());
        for (int fv = 0; fv < F[f].size(); fv++) {
            fVertList[fv] = V[F[f][fv]];
        }
        Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertList);
        // Get vector area normal
        Eigen::Vector3d fNormal = fNormals[f];

        // Build local 2D basis
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        Utils::buildPlaneBasis(fNormal, t1, t2);
        // Project p onto Newell plane and get its 2D coordinate
        Eigen::Vector3d v_proj3d = Utils::projectPointOntoPlane(fNormal, barycenter, v);
        Eigen::Vector2d v_proj2d = Utils::convertTo2D(v_proj3d, barycenter, t1, t2);

        // Project face vertices onto Newell plane using basis vectors
        std::vector<Eigen::Vector2d> fVert2D(F[f].size());
        for (int fv = 0; fv < fVertList.size(); fv++) {
            Eigen::Vector3d fv_proj3D = Utils::projectPointOntoPlane(fNormal, barycenter, fVertList[fv]);
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

        Eigen::Vector3d v_proj = Utils::revertTo3D(v_cp, barycenter, t1, t2);

        // True distance in 3D from query point
        double dist = (v - v_proj).norm();

        if (dist < min_dist) {
            min_dist = dist;
            proj = v_proj;
            elIdx = f;
        }
    }

    // Snap to nearby vertex or edge if we are too close
    // NOTE: using GLOBAL, not geodesic distance, since this is too hard for non-planar faces
    if (snap) {
        // For the face that was landed on, check if we are close to a vertex on that face
        std::vector<Eigen::Vector3d> fVertList(F[elIdx].size());
        for (int fv = 0; fv < F[elIdx].size(); fv++) {
            fVertList[fv] = V[F[elIdx][fv]];
        }
        
        // Vertex check
        for (int fv = 0; fv < F[elIdx].size(); fv++) {
            // Once found, we can just return immediately
            if ((proj - fVertList[fv]).norm() <= tol) {
                proj = fVertList[fv];
                elIdx = F[elIdx][fv];
                return 0;
            }
        }

        // If not, check if we are close to an edge in 3D
        for (int fv = 0; fv < F[elIdx].size(); fv++) {
            int fv1 = (fv + 1) % F[elIdx].size();
            Eigen::Vector3d v_projE = Utils::closestPointOnSegment3D(proj, fVertList[fv], fVertList[fv1]);
            // Once found, we can just return immediately
            if ((proj - v_projE).norm() <= tol) {
                proj = fVertList[fv];
                elIdx = heMesh.edgeIdxFromVerts(F[elIdx][fv], F[elIdx][fv1]);
                return 1;
            }
        }
    }

    // If not snapping, then we must be on a face
    // Check if we are on a non-planar face. If so, pin-point the location using MVC
    Eigen::VectorXd fHeight = H[elIdx];
    bool planar = true;
    for (int v = 0; v < fHeight.size(); v++) {
        if (std::abs(fHeight(v)) >= 1e-6) {
            planar = false;
        }
    }
    if (planar || fHeight.size() == 3) {   // Planar face, no MVC interpolation to be done
        return 2;
    }

    // Otherwise, we need to compute mean value coordinates to get projection location
    // First, project all points into Newell plane
    std::vector<Eigen::Vector3d> fVertList(F[elIdx].size());
    for (int fv = 0; fv < F[elIdx].size(); fv++) {
        fVertList[fv] = V[F[elIdx][fv]];
    }
    Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertList);
    // Get vector area normal
    Eigen::Vector3d fNormal = fNormals[elIdx];
    // Build local 2D basis
    Eigen::Vector3d t1;
    Eigen::Vector3d t2;
    Utils::buildPlaneBasis(fNormals[elIdx], t1, t2);
    Eigen::Vector2d v_proj2d = Utils::convertTo2D(proj, barycenter, t1, t2);

    std::vector<Eigen::Vector2d> fVert2D(F[elIdx].size());
    for (int fv = 0; fv < fVertList.size(); fv++) {
        Eigen::Vector3d fv_proj3D = Utils::projectPointOntoPlane(fNormal, barycenter, fVertList[fv]);
        fVert2D[fv] = Utils::convertTo2D(fv_proj3D, barycenter, t1, t2);
    }
    // Second, apply mean value coordinates to get height function weights
    Eigen::VectorXd MVCWeights(fHeight.size());
    Utils::meanValueCoordinates(v_proj2d, fVert2D, MVCWeights);

    // Now recover the height using MVC weights
    double h = MVCWeights.dot(fHeight);
    // Add height to current Newell projection
    proj += h * fNormal;

    return 2;
}

}   // namespace Mesh