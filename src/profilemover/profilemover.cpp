#include "profilemover.hpp"

#include "mesh/mesh.hpp"
#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>


namespace ProfileMover {
    profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF, 
                               const std::vector<Eigen::Vector3d> Controls, const std::vector<Eigen::Vector3d> Tangents, 
                               const std::vector<std::array<int, 4>> Splines, int alpha = 5) {
        applyMesh(meshV, meshF);
        applyCurvenet(Controls, Tangents, Splines, alpha);
        
        precomputation();
    }

    void profilemover::applyMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF) {
        if (!M_init) {
            M = Mesh::mesh(meshV, meshF);
            M_init = true;
        }
    }

    void profilemover::applyCurvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, int alpha = 5) {
        if (M_init && !CN_init) {
            CN = Curvenet::curvenet(Controls, Tangents, Splines, &M, alpha);
            CN_init = true;
        } else {
            // throw an error
        }
        return;
    }

    void profilemover::precomputation() {
        if (!M_init || !CN_init) {
            // Throw an error
            return;
        } else if (!dCN_init) {
            // Initialize discrete curvenet
            dCN = DCurvenet::dcurvenet(&CN, M.getMeanE());
            dCN_init = true;
        }
        // Compute Cut-mesh
        CM = Mesh::cutmesh(&M, &dCN);
        CM_init = true;

        // Compute operators
        CM.computeHEMap(heToCMhe, CMheTohe);
        V = CM.computeVMatrix(vToCM, mToV, heToCMhe);
        C = CM.computeCMatrix(cToCM, mToC, heToCMhe);
        // Face-based Laplacian with values pushed onto halfedges
        Eigen::SparseMatrix<double> L = CM.computeHELaplacian(CMheTohe);
        // mVtL
        mVtL = -1 * V.transpose() * L;
        // Factor V^TLV
        Eigen::SparseMatrix<double> VtLV_Mat = V.transpose() * L * V;
        VtLV.analyzePattern(VtLV_Mat);
        VtLV.factorize(VtLV_Mat);
        return;
    }

    // Runtime solvers
    // Runtime deformation
    std::vector<Eigen::Vector3d> profilemover::deform(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents) {
        if (!M_init || !CN_init || !dCN_init || !CM_init) {
            // TODO: throw error
        }
        // 1. Compute new curvenet
        CN.updateCurveNet(Controls, Tangents);
        // 2. Compute new discrete curvenet and frames
        dCN.updateDiscCurveNet();
        // FIRST SOLVE: Deformation gradients
        // Compute flattened deformation gradient matrix
        Eigen::MatrixXd f_c = CM.computeDefGrads(cToCM);    // TODO: This can be done in parallel over halfedges
        // Solve system to get interpolated 
        Eigen::MatrixXd f_v = VtLV.solve(mVtL * (C * f_c));
        // Fold back together and redistribute to their vertices
        CM.applyDefGrads(f_v, vToCM);                   // TODO: This can be done in parallel
        // SECOND SOLVE: Positions
        // Estimate new projected positions using the distributed def grads
        Eigen::MatrixXd x_c = CM.estimateCNPositions(cToCM);
        // Compute per-face deformation matrix + assemble
        Eigen::MatrixXd y_h = CM.estimateFaceDeformations(CMheTohe);
        // Compute new positions
        Eigen::MatrixXd x_v = VtLV.solve(mVtL * (C * x_c - y_h));

        return assembleFinalPositions(x_v, x_c);
    }

    // Assemble final positions into our standard data type
    // TODO: Can we not precompute a sparse operator which does this automatically, and then simply transform result to a std::vector?
    std::vector<Eigen::Vector3d> profilemover::assembleFinalPositions(Eigen::MatrixXd x_v, Eigen::MatrixXd x_c) {
        std::vector<Eigen::Vector3d> newV(M.getNumActiveV());
        // Average to get the constraint positions
        for (int m = 0; m < mToC.size(); m++) {
            Eigen::Vector3d new_v = Eigen::Vector3d::Zero();
            int m_size = mToC[m].size();
            for (int c = 0; c < m_size; c++) {
                new_v += x_c.row(mToC[m][c]).transpose();
            }
            new_v /= m_size;
            newV[m] = new_v;
        }
        // Directly copy to get the new position for original mesh vertices
        for (int v = 0; v < mToV.size(); v++) {
            newV[v] = x_v.row(mToV[v]).transpose();
        }
        return newV;
    }

}   // namespace ProfileMover