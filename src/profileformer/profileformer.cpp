#include "profileformer.hpp"

#include <limits>
#include <stdexcept>


namespace ProfileFormer {

namespace {
Eigen::Vector3d safeNormalize(const Eigen::Vector3d& v, const Eigen::Vector3d& fallback) {
    const double n = v.norm();
    if (n <= std::numeric_limits<double>::epsilon()) {
        return fallback;
    }
    return v / n;
}
}  // namespace

profileformer::profileformer(const std::vector<Eigen::Vector3d>& meshV,
                             const std::vector<std::vector<int>>& meshT,
                             const std::vector<Eigen::Vector3d>& controlP,
                             const std::vector<std::vector<int>>& curveC) {
    initializeMesh(meshV, meshT);
    initializeCurvenet(controlP, curveC);
}

void profileformer::initializeMesh(const std::vector<Eigen::Vector3d>& meshV,
                                   const std::vector<std::vector<int>>& meshT) {
    if (meshV.empty()) {
        throw std::invalid_argument("Mesh vertex list is empty");
    }
    if (meshT.empty()) {
        throw std::invalid_argument("Mesh face list is empty");
    }
    nMeshV = meshV;
    nMeshT = meshT;
}

void profileformer::initializeCurvenet(const std::vector<Eigen::Vector3d>& controlP,
                                       const std::vector<std::vector<int>>& curveC) {
    // For this initialization pass, derive normals directly from control positions.
    std::vector<Eigen::Vector3d> surfaceN(controlP.size(), Eigen::Vector3d::UnitZ());
    for (std::size_t i = 0; i < controlP.size(); ++i) {
        surfaceN[i] = safeNormalize(controlP[i], Eigen::Vector3d::UnitZ());
    }
    nCurvenet = Curvenet::curvenet(controlP, surfaceN, curveC);
}

int profileformer::precomputation() {
    // Minimal precomputation needed for rendering a stable sampled curvenet.
    // Full cutmesh/DEC precompute remains TODO.
    return 0;
}

int profileformer::deformation(const std::vector<Eigen::Vector3d>& controlP) {
    // Deformation pipeline is not needed for current visualization-only milestone.
    nCurvenet.modifyPositions(controlP);
    return 0;
}

void profileformer::computeDeformationGradients(DCurvenet::dcurvenet&, std::vector<Eigen::VectorXd>&) {
    // Legacy stub preserved: not required for curvenet visualization initialization.
}

}  // namespace ProfileFormer
