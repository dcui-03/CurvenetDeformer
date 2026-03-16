// utils.cpp

#include "utils.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <random>

namespace Utils {

// HELPERS FOR CONVERSION/COPYING

#if 0
// Legacy GLM conversion helpers — disabled to keep core free of GLM dependency.
Eigen::Vector3d glmToEigen(const glm::vec3 input) {
    Eigen::Vector3d output;
    output(0) = static_cast<double>(input.x);
    output(1) = static_cast<double>(input.y);
    output(2) = static_cast<double>(input.z);
    return output;
}

glm::vec3 eigenToGLM(const Eigen::Vector3d input) {
    glm::vec3 output;
    output.x = static_cast<float>(input(0));
    output.y = static_cast<float>(input(1));
    output.z = static_cast<float>(input(2));
    return output;
}

void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM) {
    GLM.clear();
    GLM.resize(Eig.size());
    for (int v = 0; v < Eig.size(); v++) {
        GLM[v] = eigenToGLM(Eig[v]);
    }
}

void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM) {
    Eig.clear();
    Eig.resize(GLM.size());
    for (int v = 0; v < GLM.size(); v++) {
        Eig[v] = glmToEigen(GLM[v]);
    }
}
#endif

void copyPositions(const std::vector<Eigen::Vector3d>& oldV, std::vector<Eigen::Vector3d>& newV) { newV = oldV; }

void copyConnectivity(const std::vector<std::vector<int>>& oldT, std::vector<std::vector<int>>& newT) { newT = oldT; }

Eigen::Vector3d anyUnitTangent(const Eigen::Vector3d& normal) {
    Eigen::Vector3d n = normal.normalized();
    Eigen::Vector3d candidate = std::abs(n.z()) < 0.9 ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
    Eigen::Vector3d t = n.cross(candidate);
    if (t.norm() <= std::numeric_limits<double>::epsilon()) {
        t = n.cross(Eigen::Vector3d::UnitY());
    }
    return t.normalized();
}

Eigen::Vector3d anyPerpendicularUnit(const Eigen::Vector3d& tangent) {
    Eigen::Vector3d tn = tangent.normalized();
    Eigen::Vector3d candidate = std::abs(tn.z()) < 0.9 ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d::UnitX();
    Eigen::Vector3d p = tn.cross(candidate);
    if (p.norm() <= std::numeric_limits<double>::epsilon()) {
        p = tn.cross(Eigen::Vector3d::UnitY());
    }
    return p.normalized();
}

// VECTOR/PROJECTION HELPERS

// Find basis vectors for a planar region (ex. tangent plane)
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
bool pointInPolygon2D(const Eigen::Vector2d& p, const std::vector<Eigen::Vector2d>& poly) {
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

void buildPolyscopeCurveNetwork(const Curvenet::curvenet& cn,
                                std::vector<std::array<double, 3>>& points,
                                std::vector<std::array<std::size_t, 2>>& edges,
                                std::size_t samplesPerSpline) {
    points.clear();
    edges.clear();

    for (const auto& s : cn.getSplines()) {
        const std::vector<Eigen::Vector3d> sampled = s.sampleParameterization(static_cast<int>(samplesPerSpline));
        const std::size_t base = points.size();
        for (std::size_t i = 0; i < sampled.size(); ++i) {
            const auto& p = sampled[i];
            points.push_back({p.x(), p.y(), p.z()});
            if (i > 0) {
                edges.push_back({base + i - 1, base + i});
            }
        }
    }
}

void buildPolyscopeDiscreteCurveNetwork(const DCurvenet::dcurvenet& dcn,
                                        std::vector<std::array<double, 3>>& points,
                                        std::vector<std::array<std::size_t, 2>>& edges) {
    points.clear();
    edges.clear();

    const auto& verts = dcn.verts();
    points.reserve(verts.size());
    for (const auto& v : verts) {
        const auto& p = v.position();
        points.push_back({p.x(), p.y(), p.z()});
    }

    const auto& segs = dcn.segments();
    edges.reserve(segs.size());
    for (const auto& s : segs) {
        edges.push_back({static_cast<std::size_t>(s.startDvert()), static_cast<std::size_t>(s.endDvert())});
    }
}

void buildPolyscopeControlCornerNormals(const DCurvenet::dcurvenet& dcn,
                                        std::vector<std::array<double, 3>>& origins,
                                        std::vector<std::array<double, 3>>& vectors) {
    origins.clear();
    vectors.clear();

    const auto& verts = dcn.verts();
    for (const auto& v : verts) {
        if (!v.isControl()) {
            continue;
        }

        const auto& p = v.position();
        const auto& normals = v.cornerNormals();
        for (const auto& n : normals) {
            const double len = n.norm();
            Eigen::Vector3d nu = Eigen::Vector3d::Zero();
            if (len > std::numeric_limits<double>::epsilon()) {
                nu = n / len;
            }
            origins.push_back({p.x(), p.y(), p.z()});
            vectors.push_back({nu.x(), nu.y(), nu.z()});
        }
    }
}

std::vector<std::array<double, 3>> buildControlCloud(const Curvenet::curvenet& cn) {
    std::vector<std::array<double, 3>> out;
    out.reserve(cn.controls().size());
    for (const auto& c : cn.controls()) {
        const auto& p = c.getPosition();
        out.push_back({p.x(), p.y(), p.z()});
    }
    return out;
}

}  // namespace Utils
