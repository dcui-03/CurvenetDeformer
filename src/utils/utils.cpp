#include "utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <glm/vec3.hpp>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

namespace Utils {

// HELPERS FOR CONVERSION/COPYING

// GLM::vec3 to Eigen::Vector3d converter
Eigen::Vector3d glmToEigen(const glm::vec3 input) {
    Eigen::Vector3d output;
    output(0) = static_cast<double>(input.x);
    output(1) = static_cast<double>(input.y);
    output(2) = static_cast<double>(input.z);
    return output;
}

// Eigen::Vector3d to GLM::vec3 converter
glm::vec3 eigenToGLM(const Eigen::Vector3d input) {
    glm::vec3 output;
    output.x = static_cast<float>(input(0));
    output.y = static_cast<float>(input(1));
    output.z = static_cast<float>(input(2));
    return output;
}

// Entire mesh conversion routine Eigen to GLM
void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM) {
    GLM.clear();
    GLM.resize(Eig.size());
    for (int v = 0; v < Eig.size(); v++) {
        GLM[v] = eigenToGLM(Eig[v]);
    }
    return;
}

// Entire mesh conversion routine GLM to Eigen
void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM) {
    Eig.clear();
    Eig.resize(GLM.size());
    for (int v = 0; v < GLM.size(); v++) {
        Eig[v] = glmToEigen(GLM[v]);
    }
    return;
}

// Copy positions and connectivity into a copied container
void copyPositions(const std::vector<Eigen::Vector3d>& V_old, std::vector<Eigen::Vector3d>& V_new) {
    V_new.clear();
    V_new.resize(V_old.size());
    for (int v = 0; v < V_old.size(); v++) {
        Eigen::Vector3d new_v = {V_old[v](0), V_old[v](1), V_old[v](2)};
        V_new[v] = new_v;
    }
    return;
}

void copyConnectivity(const std::vector<std::vector<int>>& T_old, std::vector<std::vector<int>>& T_new) {
    T_new.clear();
    T_new.resize(T_old.size());
    for (int f = 0; f < T_old.size(); f++) {
        std::vector<int> f_idxs;
        for (int v = 0; v < T_old[f].size(); v++) {
            f_idxs.push_back(T_old[f][v]);
        }
        T_new[f] = f_idxs;
    }
    return;
}

// SORTING
void doubleListIdxSort(std::vector<double>& ref_List, std::vector<int>& idx_List) {
    if (ref_List.size() != idx_List.size()) {
        return;
    }

    const int n = static_cast<int>(ref_List.size());

    for (int i = 0; i < n - 1; ++i) {
        bool swapped = false;

        for (int j = 0; j < n - i - 1; ++j) {
            if (ref_List[j] > ref_List[j + 1]) {
                std::swap(ref_List[j], ref_List[j + 1]);
                std::swap(idx_List[j], idx_List[j + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }
    }
}

// Insert an integer entry in a list between two specified values
bool insertIdxBetweenPair(std::vector<int>& idxList, int a, int b, int new_idx) {
    for (int i = 0; i < idxList.size(); ++i) {
        int j = (i + 1) % idxList.size();
        if (idxList[i] == a && idxList[j] == b) {
            idxList.insert(idxList.begin() + j, new_idx);
            return true;
        }
    }
    return false;
}

// Flattens an Eigen::Matrix3d into a 9x1 row vector
// NOTE: Does so column-wise!
Eigen::VectorXd flattenMatrix3d(const Eigen::Matrix3d& F) {
    return F.reshaped();
}

// Compresses a 9x1 Eigen::VectorXd into an Eigen::Matrix3d
// Assumes column-wise storage
Eigen::Matrix3d compressVector9d(const Eigen::VectorXd& f) {
    assert(f.size() == 9);
    return Eigen::Map<const Eigen::Matrix3d>(f.data());
}

// GEOMETRY HELPERS

// Find the closest point to a triangle
Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d> triVerts, const Eigen::Vector3d p) {
    const double eps = 1e-8;
    const Eigen::Vector3d& a = triVerts[0];
    const Eigen::Vector3d& b = triVerts[1];
    const Eigen::Vector3d& c = triVerts[2];

    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);

    // Vertex region outside A
    if (d1 <= 0.0 && d2 <= 0.0) {
        return a;
    }

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);

    // Vertex region outside B
    if (d3 >= 0.0 && d4 <= d3) {
        return b;
    }

    // Edge region AB
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double denom = d1 - d3;
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, a, b, true);
        }

        const double t = d1 / denom;
        return a + t * ab;
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);

    // Vertex region outside C
    if (d6 >= 0.0 && d5 <= d6) {
        return c;
    }

    // Edge region AC
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double denom = d2 - d6;
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, a, c, true);
        }
        const double t = d2 / denom;
        return a + t * ac;
    }

    // Edge region BC
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double denom = (d4 - d3) + (d5 - d6);
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, b, c, true);
        }
        const double t = (d4 - d3) / denom;
        return b + t * (c - b);
    }

    // Inside face region
    const double denom = va + vb + vc;
    if (std::abs(denom) <= eps) {
        // Degenerate triangle fallback: closest point among three edges.
        Eigen::Vector3d pab = Utils::closestPointOnSegment3D(p, a, b, true);
        Eigen::Vector3d pac = Utils::closestPointOnSegment3D(p, a, c, true);
        Eigen::Vector3d pbc = Utils::closestPointOnSegment3D(p, b, c, true);

        double dab = (p - pab).squaredNorm();
        double dac = (p - pac).squaredNorm();
        double dbc = (p - pbc).squaredNorm();

        if (dab <= dac && dab <= dbc) return pab;
        if (dac <= dab && dac <= dbc) return pac;
        return pbc;
    }

    const double invDenom = 1.0 / denom;
    const double vBary = vb * invDenom;
    const double wBary = vc * invDenom;
    return a + vBary * ab + wBary * ac;
}

// Converts a 3D direction into an angle in the tangent plane spanned by t1, t2.
// Returns false if the projected direction is degenerate.
bool directionAngleInPlane(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& target,
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double& theta
) {
    const double eps = 1e-12;

    Eigen::Vector3d n = normal;
    if (n.squaredNorm() <= eps) {
        return false;
    }
    n.normalize();

    Eigen::Vector3d d = target - origin;
    d = d - d.dot(n) * n;

    if (d.squaredNorm() <= eps) {
        return false;
    }

    const double x = d.dot(t1);
    const double y = d.dot(t2);

    theta = std::atan2(y, x);
    if (theta < 0.0) {
        theta += 2.0 * M_PI;
    }

    return true;
}

// Returns true if two angular values are effectively the same direction.
bool anglesCoincident(double a, double b, double eps) {
    double diff = std::abs(a - b);
    diff = std::min(diff, 2.0 * M_PI - diff);
    return diff <= eps;
}

// Given two unit vectors, compute the rotation from one to the other
// Rotation formulation taken from https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
Eigen::Matrix3d computeRotation(const Eigen::Vector3d& u, const Eigen::Vector3d& v) {
    double eps = 1e-6;
    // Check for degenerate vectors
    if (u.norm() <= eps || v.norm() <= eps) {
        return Eigen::Matrix3d::Zero();
    }
    // For safety, re-normalize
    Eigen::Vector3d u_norm = u.normalized();
    Eigen::Vector3d v_norm = v.normalized();
    double cosUV = u_norm.dot(v_norm);
    if (cosUV >= 1-eps) {   // Same vector
        return Eigen::Matrix3d::Identity();
    } else if (cosUV <= -1 + eps) { // Opposite vectors
        return -1 * Eigen::Matrix3d::Identity();
    }

    // Compute angle
    Eigen::Vector3d axis = u_norm.cross(v_norm);
    double sinUV = axis.norm(); // u and v are already unit
    axis.normalize();

    // Construct the rotation
    Eigen::Matrix3d rot = Eigen::Matrix3d::Zero();
    // row 0
    rot(0, 0) = axis(0) * axis(0) * (1 - cosUV) + cosUV;
    rot(0, 1) = axis(0) * axis(1) * (1 - cosUV) - axis(2)*sinUV;
    rot(0, 2) = axis(0) * axis(2) * (1 - cosUV) + axis(1)*sinUV;
    // row 1
    rot(1, 0) = axis(1) * axis(0) * (1 - cosUV) + axis(2)*sinUV;
    rot(1, 1) = axis(1) * axis(1) * (1 - cosUV) + cosUV;
    rot(1, 2) = axis(1) * axis(2) * (1 - cosUV) - axis(0)*sinUV;
    // row 2
    rot(2, 0) = axis(2) * axis(0) * (1 - cosUV) - axis(1)*sinUV;
    rot(2, 1) = axis(2) * axis(1) * (1 - cosUV) + axis(0)*sinUV;
    rot(2, 2) = axis(2) * axis(2) * (1 - cosUV) + cosUV;
    return rot;
}

// Overload
Eigen::Matrix3d computeRotation(const Eigen::Vector3d& axis, const double& theta) {
    double eps = 1e-6;
    // Check for degenerate vectors
    if (axis.norm() <= eps) {
        return Eigen::Matrix3d::Zero();
    } else if (theta <= eps) {  // No rotation
        return Eigen::Matrix3d::Identity();
    } else if (theta >= M_PI-eps && theta <= M_PI+eps) {    // 180 rotation
        return -1 * Eigen::Matrix3d::Identity();
    }
    // For safety, re-normalize
    Eigen::Vector3d a_norm = axis.normalized();
    double cosTheta = std::cos(theta);
    double sinTheta = std::sin(theta);

    // Construct the rotation
    Eigen::Matrix3d rot = Eigen::Matrix3d::Zero();
    // row 0
    rot(0, 0) = axis(0) * axis(0) * (1 - cosTheta) + cosTheta;
    rot(0, 1) = axis(0) * axis(1) * (1 - cosTheta) - axis(2)*sinTheta;
    rot(0, 2) = axis(0) * axis(2) * (1 - cosTheta) + axis(1)*sinTheta;
    // row 1
    rot(1, 0) = axis(1) * axis(0) * (1 - cosTheta) + axis(2)*sinTheta;
    rot(1, 1) = axis(1) * axis(1) * (1 - cosTheta) + cosTheta;
    rot(1, 2) = axis(1) * axis(2) * (1 - cosTheta) - axis(0)*sinTheta;
    // row 2
    rot(2, 0) = axis(2) * axis(0) * (1 - cosTheta) - axis(1)*sinTheta;
    rot(2, 1) = axis(2) * axis(1) * (1 - cosTheta) + axis(0)*sinTheta;
    rot(2, 2) = axis(2) * axis(2) * (1 - cosTheta) + cosTheta;
    return rot;
}

// Find basis vectors for a planar region (ex. tangent plane)
// Build plane basis given only n, and unitialized t1, t2
void buildPlaneBasis(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2) {
    if (std::abs(n(0)) < 0.9)
        t1 = n.cross(Eigen::Vector3d::UnitX()).normalized();
    else
        t1 = n.cross(Eigen::Vector3d::UnitY()).normalized();
    t2 = n.cross(t1); // already unit
}

// Given a point on a plane basis and the plane basis, convert to 2D planar point
Eigen::Vector2d convertTo2D(const Eigen::Vector3d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
    Eigen::Vector3d vec = p - origin;
    return Eigen::Vector2d(vec.dot(t1), vec.dot(t2));
}

// Given a 2D planar point and the plane basis, revert to its 3D counterapart
Eigen::Vector3d revertTo3D(const Eigen::Vector2d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
    return origin + p(0) * t1 + p(1) * t2;
}

// Projects a vector onto the tangent plane of a normal vector.
// If the resulting projection is near-degenerate, then return a random unit tangent vector.
double projectVectorOntoTangentPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& vec, Eigen::Vector3d& proj, double scale) {
    const double eps = 1e-12;
    // Ensure normal is unit
    Eigen::Vector3d n = normal.normalized();

    // Project onto tangent plane (i.e., subtract out projection onto normal vector)
    proj = vec - vec.dot(n) * n;
    double len = proj.norm();

    // Under degeneracy, sample a random unit vector in the tangent plane using some
    // constructed tangent plane basis
    if (len < eps) {
        // Build orthonormal tangent basis
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        buildPlaneBasis(normal, t1, t2);

        // Pick vector by sampling a random angle in [0, 2pi)
        static std::mt19937 gen(std::random_device{}());
        static std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
        // Construct proj using basis vectors
        double theta = angle_dist(gen);
        proj = std::cos(theta) * t1 + std::sin(theta) * t2;
        return -1;
    }
    // Normalize, scale, and return length
    proj /= len;
    proj *= scale;
    return len;
}

// Projects a point onto the tangent plane of a normal given a center 
Eigen::Vector3d projectPointOntoPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& center, const Eigen::Vector3d& p) {
    return p - (p - center).dot(normal) * normal;
}

// Check if a 2D point is in a 2D polygon
// To do this, we do raycasting to the segment
bool pointInPolygon2D(const Eigen::Vector2d& p, const vector2dList& poly) {
    bool inside = false;
    int n = poly.size();

    for (int v0 = 0; v0 < n; v0++) {
        int v1 = (v0 - 1) % n;
        const Eigen::Vector2d& p0 = poly[v0];
        const Eigen::Vector2d& p1 = poly[v1];
        // If 
        bool intersect = ((p0(1) > p(1)) != (p1(1) > p(1))) &&
                         (p(0) < (p1(0) - p0(0)) * (p(1) - p0(1)) / (p1(1) - p0(1)) + p0(0));

        if (intersect) {
            inside = !inside;
        }
    }
    return inside;
}

// Get the closest point on a segment in 2D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
Eigen::Vector2d closestPointOnSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, bool clip) {
    Eigen::Vector2d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16) {
        return v0;
    }
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return v0 + t * vec;
}

bool raycastToSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& direc, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1,
                        double& t, double& u, bool clip) {
    Eigen::Vector2d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16) {
        return false;
    }
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return true;
}

// Get the closest point on a segment in 3D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
Eigen::Vector3d closestPointOnSegment3D(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, bool clip) {
    Eigen::Vector3d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16)
        return v0;
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return v0 + t * vec;
}

// Computes the nearest point where a ray intersects with a face; returns false if no intersection
// Note, we build the plane to be centered at the start point for simplicity
// el_type describes if we hit an edge (1) or a vertex (2)
// t is the scale along the ray, u is the scale along the specified edge (if we intersect an edge)
// theta is the angle at the intersection (this is different for edge and vertex intersections)
// tol should be the snapping criteria
bool computeFaceIntersection(const std::vector<Eigen::Vector3d>& fVerts, const Eigen::Vector3d& fNormal,
                             const Eigen::Vector3d& start, const Eigen::Vector3d& direc,
                             double& t, int& el_type, int& local_idx, double& u, double& theta, double tol = 1e-6) {
    // build a 2D basis
    Eigen::Vector3d t1, t2;
    buildPlaneBasis(fNormal, t1, t2);
    // Project face into 2D
    vector2dList fVerts2D(fVerts.size());
    for (int fv = 0; fv < fVerts.size(); fv++) {
        Eigen::Vector3d fVertProj = projectPointOntoPlane(fNormal, start, fVerts[fv]);
        fVerts2D[fv] = convertTo2D(fVertProj, start, t1, t2);
    }
    // Project start onto 2D: We can do this if we center our Newell face at the start
    Eigen::Vector2d start2D = {0.0, 0.0};
    Eigen::Vector3d direcProj;
    double direcLen = projectVectorOntoTangentPlane(fNormal, direc, direcProj);
    Eigen::Vector2d direc2D = convertTo2D(start+direcProj, start, t1, t2);  // TODO: maybe a cleaner way to get vector into basis.

    // Iteratively raycast on all edges to find the smallest distance
    t = std::numeric_limits<double>::infinity();
    bool fHit = true;
    for (int ei = 0; ei < fVerts2D.size(); ei++) {
        double t_ei, u_ei;
        int e_next = (ei+1)%fVerts2D.size();
        bool eHit = raycastToSegment2D(start, direc2D, fVerts2D[ei], fVerts2D[e_next], t_ei, u_ei);
        // We only track this one if it's closer than our current, and a valid intersection
        if (eHit && t_ei > tol && t_ei < t) {  // TODO: this tol is not necessarily snapping criteria
            fHit = true;
            t = t_ei;
            u = u_ei;
            // Grab vertex or edge
            // This is where vertex snapping occurs
            if (u_ei <= tol) {
                u_ei = 0.0;
                el_type = 2;
                local_idx = ei;
            } else if (u_ei >= 1 - tol) {
                u_ei = 1.0;
                el_type = 2;
                local_idx = e_next;
            } else {
                el_type = 1;
                local_idx = ei;
            }
        }
    }
    // If no hits, then just return false
    if (!fHit) {
        return false;
    }

    // Otherwise, find the angle between the direc and the CCW face edge
    Eigen::Vector2d edgeVec;
    if (el_type == 1) { // edge case
        theta = vectorAngle(start2D, start2D+direc2D, fVerts2D[local_idx], fVerts2D[(local_idx+1)%fVerts2D.size()]);
    } else {    // vertex case
        theta = vectorAngle(start2D, start2D+direc2D, fVerts2D[(local_idx + fVerts2D.size() - 1)%fVerts2D.size()], fVerts2D[local_idx]);
    }

    return true;
}

// Compute face intersection w.r.t. a target
// NOTE: Unlike the regular face intersection, returns True IFF the target is closest, and False otherwise.
// TODO: combine the two face intersection methods?
// TODO: If intersects with a vertex, compute and store the ccw edge idx. Also lift direction to the plane formed by
// intersection corner in 3D
bool computeFaceIntersectionTarget(const std::vector<Eigen::Vector3d>& fVerts, const Eigen::Vector3d& fNormal,
                             const Eigen::Vector3d& start, const Eigen::Vector3d& target,
                             double& t, int& el_type, int& local_idx, double& u, double& theta, double tol = 1e-6) {
    // build a 2D basis
    Eigen::Vector3d direc = target - start;
    Eigen::Vector3d t1, t2;
    buildPlaneBasis(fNormal, t1, t2);
    // Project face into 2D
    vector2dList fVerts2D(fVerts.size());
    for (int fv = 0; fv < fVerts.size(); fv++) {
        Eigen::Vector3d fVertProj = projectPointOntoPlane(fNormal, start, fVerts[fv]);
        fVerts2D[fv] = convertTo2D(fVertProj, start, t1, t2);
    }
    // Project start onto 2D: We can do this if we center our Newell face at the start
    Eigen::Vector2d start2D = {0.0, 0.0};
    Eigen::Vector3d direcProj;
    double direcLen = projectVectorOntoTangentPlane(fNormal, direc, direcProj);
    Eigen::Vector2d direc2D = convertTo2D(start+direcProj, start, t1, t2);  // TODO: maybe a cleaner way to get vector into basis.

    // Iteratively raycast on all edges to find the smallest distance
    t = std::numeric_limits<double>::infinity();
    bool fHit = true;
    for (int ei = 0; ei < fVerts2D.size(); ei++) {
        double t_ei, u_ei;
        int e_next = (ei+1)%fVerts2D.size();
        bool eHit = raycastToSegment2D(start, direc2D, fVerts2D[ei], fVerts2D[e_next], t_ei, u_ei);
        // We only track this one if it's closer than our current, and a valid intersection
        if (eHit && t_ei > tol && t_ei < t) {  // TODO: this tol is not necessarily snapping criteria
            fHit = true;
            t = t_ei;
            u = u_ei;
            // Grab vertex or edge
            // This is where vertex snapping occurs
            if (u_ei <= tol) {
                el_type = 2;
                local_idx = ei;
            } else if (u_ei >= 1 - tol) {
                el_type = 2;
                local_idx = e_next;
            } else {
                el_type = 1;
                local_idx = ei;
            }
        }
    }

    // Test if our target is closer than the intersection
    if (direcLen <= t + tol) {
        return true;
    }

    // Otherwise, find the angle between the direc and the CCW face edge
    Eigen::Vector2d edgeVec;
    if (el_type == 1) { // edge case
        // We can just use the Newell plane to get the angle
        theta = vectorAngle(start2D, start2D+direc2D, fVerts2D[local_idx], fVerts2D[(local_idx+1)%fVerts2D.size()]);
    } else {    // vertex case
        // We need to project back onto the original corner to get the "3D" corner angle
        Eigen::Vector3d edge_vec_p1 = (fVerts[(local_idx + 1)%fVerts.size()] - fVerts[local_idx]).normalized();
        Eigen::Vector3d edge_vec_m1 = (fVerts[(local_idx + fVerts.size() - 1)%fVerts.size()] - fVerts[local_idx]).normalized();
        Eigen::Vector3d corner_normal = (edge_vec_p1).cross(edge_vec_m1);
        Eigen::Vector3d corner_direc;
        projectVectorOntoTangentPlane(corner_normal, -1 * direc, corner_direc);
        corner_direc.normalize();
        theta = std::acos(std::clamp(corner_direc.dot(edge_vec_m1), -1.0, 1.0));
    }

    return false;
}

// MEAN VALUE COORDINATES
// Compute angle between any two 2D vectors given the four endpoints
// Vectors are computed as p1 - p0, p3 - p2
double vectorAngle(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3) {
    Eigen::Vector2d v0 = p1 - p0;
    Eigen::Vector2d v1 = p3 - p2;
    return std::atan2(v0(0)*v1(1) - v0(1)*v1(0), v0(0)*v1(0) + v0(1)*v1(1));
}

// Compute the sign of a double value
double computeSign(const double& value) {
    if (value > 0.0) {
        return 1.0;
    } else if (value < 0.0) {
        return -1.0;
    }
    return 0.0;
}

void meanValueCoordinates(const Eigen::Vector2d& target, const vector2dList& cage, Eigen::VectorXd& weights) {
    weights.setZero();

    double W = 0.0;
    std::vector<double> beta(cage.size());
    std::vector<double> gamma(cage.size());
    std::vector<double> s(cage.size());
    std::vector<double> r(cage.size());
    for (int v = 0; v < cage.size(); v++) {
        beta[v] = vectorAngle(cage[v], cage[(v+1)%cage.size()], target, cage[v]);
        gamma[v] = vectorAngle(cage[(v+1)%cage.size()], cage[v], cage[(v+1)%cage.size()], target);
        s[v] = beta[v] + gamma[v];
        r[v] = (cage[v] - target).norm();
    }
    std::vector<double> w(cage.size());
    for (int v = 0; v < cage.size(); v++) {
        int v_m1 = (v + cage.size() - 1)%cage.size();
        double alpha_m1p1 = vectorAngle(target, cage[v_m1], target, cage[(v+1)%cage.size()]);
        double s_m1p1 = M_PI * (computeSign(s[v_m1]) + computeSign(s[v])) - s[v_m1] - s[v];
        if (computeSign(alpha_m1p1) != computeSign(s_m1p1)) {    // NOTE: Potential numerical issue here
            alpha_m1p1 *= -1.0;
        }
        w[v] = r[v_m1] * std::sin(alpha_m1p1/2.0);
        for (int u = 0; u < cage.size(); u++) {
            if ((u != v_m1) && (u != v)) {
                w[v] *= r[u] * std::sin(std::abs(s[v])/2.0);
            }
        }
        W += w[v];
    }
    // Normalize to get final weights
    for (int v = 0; v < cage.size(); v++) {
        weights[v] = w[v]/W;
    }

    return;
}

} // namespace Utils