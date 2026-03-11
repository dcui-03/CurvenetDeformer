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

namespace Utils {

#if 0
// Legacy conversion helpers preserved from original stub.
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
    for (int v = 0; v < static_cast<int>(Eig.size()); v++) {
        GLM[v] = eigenToGLM(Eig[v]);
    }
}

void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM) {
    Eig.clear();
    Eig.resize(GLM.size());
    for (int v = 0; v < static_cast<int>(GLM.size()); v++) {
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

Eigen::Vector3d projectAndNormalizeToTangentPlane(const Eigen::Vector3d& normal,
                                                  const Eigen::Vector3d& tangent) {
    const Eigen::Vector3d tn = tangent.normalized();
    const Eigen::Vector3d projected = normal - tn * normal.dot(tn);
    return projected.normalized();
}

namespace {
std::uint64_t edgeKey(int a, int b) {
    const std::uint32_t lo = static_cast<std::uint32_t>(std::min(a, b));
    const std::uint32_t hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | static_cast<std::uint64_t>(hi);
}
}  // namespace

double computeMeanMeshEdgeLength(const std::vector<Eigen::Vector3d>& verts,
                                 const std::vector<std::vector<int>>& faces) {
    std::unordered_set<std::uint64_t> seen;
    double sum = 0.0;
    std::size_t count = 0;

    for (const auto& f : faces) {
        if (f.size() < 2) {
            continue;
        }
        const std::size_t m = f.size();
        for (std::size_t i = 0; i < m; ++i) {
            const int a = f[i];
            const int b = f[(i + 1) % m];
            if (a < 0 || b < 0 || static_cast<std::size_t>(a) >= verts.size() || static_cast<std::size_t>(b) >= verts.size() ||
                a == b) {
                continue;
            }

            const std::uint64_t key = edgeKey(a, b);
            if (!seen.insert(key).second) {
                continue;
            }

            sum += (verts[static_cast<std::size_t>(a)] - verts[static_cast<std::size_t>(b)]).norm();
            ++count;
        }
    }

    if (count == 0) {
        return 1.0;
    }
    return sum / static_cast<double>(count);
}

void loadObjMesh(const std::string& path,
                 std::vector<Eigen::Vector3d>& vertices,
                 std::vector<std::vector<int>>& faces) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open OBJ: " + path);
    }

    vertices.clear();
    faces.clear();

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string tag;
        iss >> tag;

        if (tag == "v") {
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            iss >> x >> y >> z;
            vertices.emplace_back(x, y, z);
        } else if (tag == "f") {
            std::vector<int> face;
            std::string tok;
            while (iss >> tok) {
                const std::size_t slash = tok.find('/');
                const std::string idxStr = (slash == std::string::npos) ? tok : tok.substr(0, slash);
                if (idxStr.empty()) {
                    continue;
                }

                const int idxRaw = std::stoi(idxStr);
                int idx = 0;
                if (idxRaw > 0) {
                    idx = idxRaw - 1;
                } else {
                    const int fromBack = -idxRaw;
                    if (fromBack <= 0 || static_cast<std::size_t>(fromBack) > vertices.size()) {
                        throw std::runtime_error("Invalid negative face index in OBJ");
                    }
                    idx = static_cast<int>(vertices.size()) - fromBack;
                }
                face.push_back(idx);
            }

            if (face.size() >= 3) {
                faces.push_back(std::move(face));
            }
        }
    }
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
