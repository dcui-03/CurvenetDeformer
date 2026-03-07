// utils.cpp

#include "utils.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
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
