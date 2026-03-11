// profileformer.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#if 0
// Legacy includes preserved from original stub.
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#endif

namespace ProfileFormer {

class profileformer {
public:
    // Constructor, which takes the mesh vertices and connectivity, as well as the spline controls and their connectivities
    // NOTE: We define each curve by [start, tangent1, tangent2, end]
    profileformer(const std::vector<Eigen::Vector3d>& meshV,
                  const std::vector<std::vector<int>>& meshT,
                  const std::vector<Eigen::Vector3d>& controlP,
                  const std::vector<Eigen::Vector3d>& surfaceN,
                  const std::vector<std::vector<int>>& curveC);

    // Precompute cut-mesh and operators
    int precomputation();

    // Apply deformation given the new control point locations (connectivity should be same)
    int deformation(const std::vector<Eigen::Vector3d>& controlP);

    // TODO: Initialize directly from existing mesh data struct (ex. GeometryCentral or minimesh)
    // TODO: Needs getters so that the UI can ask for internals
    const Curvenet::curvenet& getNeutralCurvenet() const { return nCurvenet; }
    const DCurvenet::dcurvenet& getNeutralDCurvenet() const { return nDCurvenet; }
    const std::vector<Eigen::Vector3d>& getNeutralMeshV() const { return nMeshV; }
    const std::vector<std::vector<int>>& getNeutralMeshT() const { return nMeshT; }

protected:
    // No class inheritance
private:
    // initialization functions called by constructor
    // TODO: We should probably compute vertex, edge, and face normals here and store them somewhere...
    void initializeMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshT);
    void initializeCurvenet(const std::vector<Eigen::Vector3d>& controlP,
                            const std::vector<Eigen::Vector3d>& surfaceN,
                            const std::vector<std::vector<int>>& curveC);
    void initializeDCurvenet();

    // Computes deformations on the temporary discrete curvenet, then fills in a list of flattened deformation
    // gradients per vertex
    void computeDeformationGradients(DCurvenet::dcurvenet& tempDC, std::vector<Eigen::VectorXd>& defGrads);

    // Store copy of initial mesh
    // TODO: Need to pick a HE mesh class (ex. GeometryCentral or minimesh)
    std::vector<Eigen::Vector3d> nMeshV;
    std::vector<std::vector<int>> nMeshT;

    // NOTE: no need to store updated states curvenet and dcurvenet, we have to create new copies at execution time
    // Store the neutral curvenet (spline)
    Curvenet::curvenet nCurvenet;
    // Store the neutral discrete curvenet
    DCurvenet::dcurvenet nDCurvenet;

    // Discretization parameters / cached neutral statistics.
    int dcurve_samples_per_mean_edge = 5;
    int dcurve_uniform_refine_samples = 64;
    double neutral_mean_edge_length = 1.0;

    // Store operators (left as placeholders; not required for curvenet visualization init)
    Eigen::SparseMatrix<double> VtLV;
    Eigen::SparseMatrix<double> mVtL;

    Eigen::SparseMatrix<double> L;
    Eigen::SparseMatrix<double> V;
    Eigen::SparseMatrix<double> C;
};

}  // namespace ProfileFormer
