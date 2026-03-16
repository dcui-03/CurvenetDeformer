// profilemover.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "dcurvenet/pdcurvenet.hpp"
#include "mesh/mesh.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <memory>
#include <vector>

namespace ProfileMover {

class profilemover {
public:
    profilemover() = default;
    profilemover(const std::vector<Eigen::Vector3d>& meshV,
                 const std::vector<std::vector<int>>& meshT,
                 const std::vector<Eigen::Vector3d>& controlP,
                 const std::vector<Eigen::Vector3d>& surfaceN,
                 const std::vector<std::vector<int>>& curveC);
    profilemover(const std::vector<Eigen::Vector3d>& meshV,
                 const std::vector<std::vector<int>>& meshT,
                 const Curvenet::curvenet& neutralCurvenet);

    int precomputation();
    int deformation(const std::vector<Eigen::Vector3d>& controlP);

    const Curvenet::curvenet& getNeutralCurvenet() const { return nCurvenet; }
    const DCurvenet::dcurvenet& getNeutralDCurvenet() const { return nDCurvenet; }
    const DCurvenet::pdcurvenet& getNeutralPDCurvenet() const { return nPDCurvenet; }
    const std::vector<Eigen::Vector3d>& getNeutralMeshV() const { return nMeshV; }
    const std::vector<std::vector<int>>& getNeutralMeshT() const { return nMeshT; }
    Mesh::mesh& getNeutralMesh() { return *nMesh; }
    int getDCurveSamplesPerMeanEdge() const { return dcurve_samples_per_mean_edge; }
    int getDCurveUniformRefineSamples() const { return dcurve_uniform_refine_samples; }
    void setDiscretizationParameters(int samplesPerMeanEdge, int uniformRefineSamples);

private:
    void initializeMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshT);
    void initializeCurvenet(const std::vector<Eigen::Vector3d>& controlP,
                            const std::vector<Eigen::Vector3d>& surfaceN,
                            const std::vector<std::vector<int>>& curveC);
    void initializeCurvenet(const Curvenet::curvenet& neutralCurvenet);
    void initializeDCurvenet();
    void initializePDCurvenet();
    void computeDeformationGradients(DCurvenet::dcurvenet& tempDC, std::vector<Eigen::VectorXd>& defGrads);

    std::vector<Eigen::Vector3d> nMeshV;
    std::vector<std::vector<int>> nMeshT;
    std::unique_ptr<Mesh::mesh> nMesh;

    Curvenet::curvenet nCurvenet;
    DCurvenet::dcurvenet nDCurvenet;
    DCurvenet::pdcurvenet nPDCurvenet;

    int dcurve_samples_per_mean_edge = 5;
    int dcurve_uniform_refine_samples = 64;
    double neutral_mean_edge_length = 1.0;

    Eigen::SparseMatrix<double> VtLV;
    Eigen::SparseMatrix<double> mVtL;
    Eigen::SparseMatrix<double> L;
    Eigen::SparseMatrix<double> V;
    Eigen::SparseMatrix<double> C;
};

}  // namespace ProfileMover
