#include "profilemover.hpp"

#include "mesh/mesh.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>

namespace ProfileMover {

profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV,
                           const std::vector<std::vector<int>>& meshT,
                           const std::vector<Eigen::Vector3d>& controlP,
                           const std::vector<Eigen::Vector3d>& surfaceN,
                           const std::vector<std::vector<int>>& curveC) {
    initializeMesh(meshV, meshT);
    initializeCurvenet(controlP, surfaceN, curveC);
    initializeDCurvenet();
    initializePDCurvenet();
}

profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV,
                           const std::vector<std::vector<int>>& meshT,
                           const Curvenet::curvenet& neutralCurvenet) {
    initializeMesh(meshV, meshT);
    initializeCurvenet(neutralCurvenet);
    initializeDCurvenet();
    initializePDCurvenet();
}

void profilemover::initializeMesh(const std::vector<Eigen::Vector3d>& meshV,
                                  const std::vector<std::vector<int>>& meshT) {
    if (meshV.empty()) {
        throw std::invalid_argument("Mesh vertex list is empty");
    }
    if (meshT.empty()) {
        throw std::invalid_argument("Mesh face list is empty");
    }
    nMeshV = meshV;
    nMeshT = meshT;
    nMesh = std::make_unique<Mesh::mesh>(nMeshV, nMeshT);
    neutral_mean_edge_length = nMesh->getMeanE();
}

void profilemover::initializeCurvenet(const std::vector<Eigen::Vector3d>& controlP,
                                      const std::vector<Eigen::Vector3d>& surfaceN,
                                      const std::vector<std::vector<int>>& curveC) {
    std::vector<Eigen::Vector3d> normals(controlP.size(), Eigen::Vector3d::UnitZ());
    for (std::size_t i = 0; i < controlP.size(); ++i) {
        normals[i] = (i < surfaceN.size() ? surfaceN[i] : controlP[i]).normalized();
    }
    nCurvenet = Curvenet::curvenet(controlP, normals, curveC);
}

void profilemover::initializeCurvenet(const Curvenet::curvenet& neutralCurvenet) {
    nCurvenet = neutralCurvenet;
}

void profilemover::initializeDCurvenet() {
    nDCurvenet = DCurvenet::dcurvenet(
        nCurvenet,
        neutral_mean_edge_length,
        dcurve_samples_per_mean_edge,
        dcurve_uniform_refine_samples);
}

void profilemover::initializePDCurvenet() {
    nPDCurvenet = DCurvenet::pdcurvenet(nDCurvenet, nMeshV, nMeshT);
}

void profilemover::setDiscretizationParameters(int samplesPerMeanEdge, int uniformRefineSamples) {
    dcurve_samples_per_mean_edge = std::max(1, samplesPerMeanEdge);
    dcurve_uniform_refine_samples = std::max(8, uniformRefineSamples);
}

int profilemover::precomputation() {
    initializeDCurvenet();
    initializePDCurvenet();
    return 0;
}

int profilemover::deformation(const std::vector<Eigen::Vector3d>& controlP) {
    nCurvenet.modifyPositions(controlP);
    return 0;
}

void profilemover::computeDeformationGradients(DCurvenet::dcurvenet&, std::vector<Eigen::VectorXd>&) {
    // Deformation pipeline not implemented yet.
}

}  // namespace ProfileMover
