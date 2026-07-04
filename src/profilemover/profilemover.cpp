#include "profilemover.hpp"

#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>


namespace ProfileMover {

    profilemover::profilemover(std::vector<Eigen::Vector3d>& meshV, std::vector<std::vector<int>>& meshF) {
        M = Mesh::mesh(meshV, meshF);
    }

    void profilemover::precomputation(std::vector<Eigen::Vector3d> Controls, 
                                      std::vector<Eigen::Vector3d> Tangents, 
                                      std::vector<std::array<int, 4>> Splines,
                                      int alpha) {
        // Initialize curvenet
        CN = Curvenet::curvenet(Controls, Tangents, Splines, M, alpha);
        // Initialize discrete curvenet
        dCN = DCurvenet::dcurvenet(&CN, M.getMeanE(), alpha);
        dCN_init = true;
        // Compute Cut-mesh
        M.computeCutMesh();
        // Compute operators
        // C
        // V
        std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> VC = M.computeVC();
        V = VC.first;
        C = VC.second;
        // Face-based Laplacian with values pushed onto halfedges
        Eigen::SparseMatrix<double> L = M.computeHELaplacian();
        // mVtL
        mVtL = -1 * V.transpose() * L;
        // Factor V^TLV
        Eigen::SparseMatrix<double> VtLV_Mat = V.transpose() * L * V;
        VtLV.analyzePattern(VtLV_Mat);
        VtLV.factorize(VtLV_Mat);
    }

    // Runtime deformation
    std::vector<Eigen::Vector3d> profilemover::deform(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents) {
        // 1. Compute new curvenet
        CN.updateCurveNet(Controls, Tangents);
        // 2. Compute new discrete curvenet and frames
        dCN.updateDiscCurveNet();
        // FIRST SOLVE: Deformation gradients
        // 3. Compute flattened deformation gradient matrix
        Eigen::MatrixXd cnDefGrads = M.computeDefGrads();
        // 4. Solve system using precomputed factorization
        // First, build RHS (TODO: This is a placeholder)
        Eigen::MatrixXd RHSdefGrad = mVtL * (cnDefGrads);
        Eigen::MatrixXd defGrads = VtLV.solve(RHSdefGrad);
        // SECOND SOLVE: Positions
        // 5. Redistribute def grads onto mesh vertices

        // 6. Compute per-face matrix + assemble

        // 7. Compute vertex projections
        Eigen::MatrixXd projDefs = M.estimateProjectionDefs();
        // 8. Assemble RHS

        // 9. Compute new positions
        Eigen::MatrixXd newPositions = VtLV.solve(projDefs);
        // Re-format to return type (maybe just return the new positions and handle at the hand-off?)
        std::vector<Eigen::Vector3d> 
        for (int i = 0; i < num_v; i++) {

        }
    }


}   // namespace ProfileMover