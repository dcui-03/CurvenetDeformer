#include "profileformer.hpp"
#include "utils/utils.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ProfileFormer {


profileformer::profileformer(const std::vector<Eigen::Vector3d>& meshV,
                             const std::vector<std::vector<int>>& meshT,
                             const std::vector<Eigen::Vector3d>& controlP,
                             const std::vector<Eigen::Vector3d>& surfaceN,
                             const std::vector<std::vector<int>>& curveC) {
    initializeMesh(meshV, meshT);
    initializeCurvenet(controlP, surfaceN, curveC);
    initializeDCurvenet();
    initializePDCurvenet();
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
    neutral_mean_edge_length = Utils::computeMeanMeshEdgeLength(nMeshV, nMeshT);
}

void profileformer::initializeCurvenet(const std::vector<Eigen::Vector3d>& controlP,
                                       const std::vector<Eigen::Vector3d>& surfaceN,
                                       const std::vector<std::vector<int>>& curveC) {
    std::vector<Eigen::Vector3d> normals(controlP.size(), Eigen::Vector3d::UnitZ());
    for (std::size_t i = 0; i < controlP.size(); ++i) {
        normals[i] = (i < surfaceN.size() ? surfaceN[i] : controlP[i]).normalized();
    }
    nCurvenet = Curvenet::curvenet(controlP, normals, curveC);
}

void profileformer::initializeDCurvenet() {
    nDCurvenet = DCurvenet::dcurvenet(
        nCurvenet,
        neutral_mean_edge_length,
        dcurve_samples_per_mean_edge,
        dcurve_uniform_refine_samples);
}

void profileformer::initializePDCurvenet() {
    nPDCurvenet = DCurvenet::pdcurvenet(nDCurvenet, nMeshV, nMeshT);
}

void profileformer::setDiscretizationParameters(int samplesPerMeanEdge, int uniformRefineSamples) {
    dcurve_samples_per_mean_edge = std::max(1, samplesPerMeanEdge);
    dcurve_uniform_refine_samples = std::max(8, uniformRefineSamples);
}

int profileformer::precomputation() {
    // Keep this cheap for UI refreshes.
    initializeDCurvenet();
    initializePDCurvenet();
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
