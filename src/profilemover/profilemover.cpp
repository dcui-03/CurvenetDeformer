#include "profilemover.hpp"

#include "mesh/mesh.hpp"
#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>


namespace ProfileMover {
    profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF, 
                               const std::vector<Eigen::Vector3d> Controls, const std::vector<Eigen::Vector3d> Tangents, 
                               const std::vector<std::array<int, 4>> Splines, int alpha) {
        applyMesh(meshV, meshF);
        std::cout << "Mesh created." << std::endl;
        applyCurvenet(Controls, Tangents, Splines, alpha);
        std::cout << "Curvenet created." << std::endl;
        computeDiscreteCurvenet();
        std::cout << "Discrete Curvenet created." << std::endl;
        computeCutMesh();
        std::cout << "Cutmesh created." << std::endl;
        
        precomputation();
        std::cout << "Operators computed." << std::endl;
    }
    // Only apply mesh
    profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF) {
        applyMesh(meshV, meshF);
    }
    // Blank init
    profilemover::profilemover() {

    }

    void profilemover::applyMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF) {
        if (M_init) {
            throw std::runtime_error("profilemover::applyMesh(): mesh already initialized");
        }
        M = Mesh::mesh(meshV, meshF);
        M_init = true;
        return;
    }

    void profilemover::applyCurvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, int alpha) {
        if (!M_init) {
            throw std::runtime_error("profilemover::applyCurvenet(): mesh must be initialized first");
            return;
        } else if (CN_init) {
            throw std::runtime_error("profilemover::applyCurvenet(): curvenet already initialized");
            return;
        }
        CN = Curvenet::curvenet(Controls, Tangents, Splines, M, alpha);
        CN_init = true;
        return;
    }

    void profilemover::computeDiscreteCurvenet() {
        if (!M_init || !CN_init) {
            throw std::runtime_error("profilemover::computeDiscreteCurvenet(): curvenet or mesh not initialized");
            return;
        } else if (dCN_init) {
            throw std::runtime_error("profilemover::computeDiscreteCurvenet(): discrete curvenet already computed");
            return;
        } else {
            // Initialize discrete curvenet
            dCN = DCurvenet::dcurvenet(&CN);
            dCN_init = true;
        }
        return;
    }

    void profilemover::computeCutMesh() {
        if (!M_init || !CN_init || !dCN_init) {
            throw std::runtime_error("profilemover::computeCutMesh(): curvenet, mesh, or discrete curvenet not initialized");
            return;
        } else if (CM_init) {
            throw std::runtime_error("profilemover::computeCutMesh(): cutmesh already computed");
            return;
        } else {
            // Compute Cut-mesh
            CM = Mesh::cutmesh(&M, &dCN);
            CM_init = true;
        }
        return;
    }

    void profilemover::precomputation() {
        if (!M_init || !CN_init || !dCN_init || !CM_init) {
            throw std::runtime_error("profilemover::precomputation(): curvenet, mesh, cutmesh, or discrete curvenet not initialized.");
            return;
        }

        // Compute operators
        CM.computeHEMap(heToCMhe, CMheTohe);
        int V_success = CM.computeVMatrix(V, vToCM, mToV, heToCMhe);
        if (V_success != 1) {
            throw std::runtime_error("profilemover::precomputation(): Unable to construct matrix V.");
            return;
        }
        int C_success = CM.computeCMatrix(C, cToCM, mToC, heToCMhe);
        if (C_success != 1) {
            throw std::runtime_error("profilemover::precomputation(): Unable to construct matrix C.");
            return;
        }
        // Face-based Laplacian with values pushed onto halfedges
        Eigen::SparseMatrix<double> L;
        int L_success = CM.computeHELaplacian(L, CMheTohe);
        if (L_success != 1) {
            throw std::runtime_error("profilemover::precomputation(): Error constructing halfedge Laplacian.");
            return;
        }
        // mVtL
        mVtL = -1 * V.transpose() * L;
        // Factor V^TLV
        Eigen::SparseMatrix<double> VtLV_Mat = V.transpose() * L * V;
        VtLV.analyzePattern(VtLV_Mat);
        VtLV.factorize(VtLV_Mat);
        return;
    }

    const Mesh::mesh& profilemover::mesh() const {
        if (!M_init) {
            throw std::runtime_error("profilemover::mesh(): mesh not initialized");
        }
        return M;
    }

    const Mesh::cutmesh& profilemover::cutmesh() const {
        if (!CM_init) {
            throw std::runtime_error("profilemover::cutmesh(): cutmesh not initialized");
        }
        return CM;
    }

    const Curvenet::curvenet& profilemover::curvenet() const {
        if (!CN_init) {
            throw std::runtime_error("profilemover::curvenet(): curvenet not initialized");
        }
        return CN;
    }

    const DCurvenet::dcurvenet& profilemover::discreteCurvenet() const {
        if (!dCN_init) {
            throw std::runtime_error("profilemover::discreteCurvenet(): discrete curvenet not initialized");
        }
        return dCN;
    }

    // Runtime solvers
    // Runtime deformation
    std::vector<Eigen::Vector3d> profilemover::deform(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents) {
        std::vector<Eigen::Vector3d> temp;
        if (!M_init || !CN_init || !dCN_init || !CM_init) {
            throw std::runtime_error("profilemover::deform(): mesh, curvenet, discrete curvenet, or cutmesh not computed");
            return temp;
        }
        // 1. Compute new curvenet
        CN.updateCurveNet(Controls, Tangents);
        // 2. Compute new discrete curvenet and frames
        dCN.updateDiscCurveNet();
        // FIRST SOLVE: Deformation gradients
        // Compute flattened deformation gradient matrix
        Eigen::MatrixXd f_c;
        int DG_success = CM.computeDefGrads(f_c, cToCM);    // TODO: This can be done in parallel over halfedges
        if (DG_success != 1) {
            throw std::runtime_error("profilemover::deform(): deformation gradient computation failed.");
            return temp;
        }
        // Solve system to get interpolated 
        Eigen::MatrixXd f_v = VtLV.solve(mVtL * (C * f_c));
        // Fold back together and redistribute to their vertices
        CM.applyDefGrads(f_v, vToCM);                   // TODO: This can be done in parallel
        // SECOND SOLVE: Positions
        // Estimate new projected positions using the distributed def grads
        Eigen::MatrixXd x_c;
        int cnPos_success = CM.estimateCNPositions(x_c, cToCM);
        if (cnPos_success != 1) {
            throw std::runtime_error("profilemover::deform(): curve network projection estimation failed.");
            return temp;
        }
        // Compute per-face deformation matrix + assemble
        Eigen::MatrixXd y_h;
        int faceDef_success = CM.estimateFaceDeformations(y_h, CMheTohe);
        if (faceDef_success != 1) {
            throw std::runtime_error("profilemover::deform(): face deformation estimate failed.");
            return temp;
        }
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