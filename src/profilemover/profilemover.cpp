#include "profilemover.hpp"

#include "mesh/mesh.hpp"
#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <chrono>
#include <iostream>


namespace ProfileMover {
    profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF, 
                               const std::vector<Eigen::Vector3d> Controls, const std::vector<Eigen::Vector3d> Tangents, 
                               const std::vector<std::array<int, 4>> Splines, int alpha, bool arap): arap(arap) {
        applyMesh(meshV, meshF);

        applyCurvenet(Controls, Tangents, Splines, alpha);
        
        computeDiscreteCurvenet();

        computeCutMesh();
        
        precomputation();
    }
    // Only apply mesh
    profilemover::profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF) {
        applyMesh(meshV, meshF);
    }
    // Blank init
    profilemover::profilemover() {

    }

    void profilemover::toggleARAP(bool toggle) {
        arap = toggle;
        return;
    }

    void profilemover::applyMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF) {
        if (M_init) {
            throw std::runtime_error("profilemover::applyMesh(): mesh already initialized");
        }
        auto mesh_start = std::chrono::steady_clock::now();
        M = Mesh::mesh(meshV, meshF);
        M_init = true;
        auto mesh_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> mesh_elapsed = mesh_end - mesh_start;
        std::cout << "Mesh created in time " << mesh_elapsed.count() << " seconds"<< std::endl;
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
        auto CN_start = std::chrono::steady_clock::now();
        CN = Curvenet::curvenet(Controls, Tangents, Splines, M, alpha);
        CN_init = true;
        auto CN_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> CN_elapsed = CN_end - CN_start;
        std::cout << "Curvenet created in time " << CN_elapsed.count() << " seconds" << std::endl;
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
            auto dCN_start = std::chrono::steady_clock::now();
            // Initialize discrete curvenet
            dCN = DCurvenet::dcurvenet(&CN);
            dCN_init = true;
            auto dCN_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> dCN_elapsed = dCN_end - dCN_start;
            std::cout << "Discrete Curvenet created in time " << dCN_elapsed.count() << " seconds" << std::endl;
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
            auto CM_start = std::chrono::steady_clock::now();
            CM = Mesh::cutmesh(&M, &dCN);
            CM_init = true;
            auto CM_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> CM_elapsed = CM_end - CM_start;
            std::cout << "Cutmesh created in time " << CM_elapsed.count() << " seconds" << std::endl;
        }
        return;
    }

    void profilemover::precomputation() {
        if (!M_init || !CN_init || !dCN_init || !CM_init) {
            throw std::runtime_error("profilemover::precomputation(): curvenet, mesh, cutmesh, or discrete curvenet not initialized.");
            return;
        }
        auto Op_start = std::chrono::steady_clock::now();
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
        auto Op_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> Op_elapsed = Op_end - Op_start;
        std::cout << "Operators computed in time " << Op_elapsed.count() << " seconds" << std::endl;
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
        auto start = std::chrono::steady_clock::now();
        CN.updateCurveNet(Controls, Tangents);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        // std::cout << "Updated Curvenet in " << elapsed.count() << " seconds"<< std::endl;
        // 2. Compute new discrete curvenet and frames
        start = std::chrono::steady_clock::now();
        dCN.updateDiscCurveNet();
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Updated discrete Curvenet in " << elapsed.count() << " seconds"<< std::endl;
        // FIRST SOLVE: Deformation gradients
        // Compute flattened deformation gradient matrix
        start = std::chrono::steady_clock::now();
        Eigen::MatrixXd f_c;
        int DG_success = CM.computeDefGrads(f_c, cToCM);    // TODO: This can be done in parallel over halfedges
        if (DG_success != 1) {
            throw std::runtime_error("profilemover::deform(): deformation gradient computation failed.");
            return temp;
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Computed cut-vert def grads in " << elapsed.count() << " seconds"<< std::endl;

        // Solve system to get interpolated def grads
        start = std::chrono::steady_clock::now();
        Eigen::MatrixXd f_v = VtLV.solve(mVtL * (C * f_c));
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Computed all def grads in " << elapsed.count() << " seconds"<< std::endl;

        // Fold back together and redistribute to their vertices
        start = std::chrono::steady_clock::now();
        CM.applyDefGrads(f_v, vToCM);                   // TODO: This can be done in parallel
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Applied all def grads in " << elapsed.count() << " seconds"<< std::endl;

        // SECOND SOLVE: Positions
        // Estimate new projected positions using the distributed def grads
        start = std::chrono::steady_clock::now();
        Eigen::MatrixXd x_c;
        int cnPos_success = CM.estimateCNPositions(x_c, cToCM);
        if (cnPos_success != 1) {
            throw std::runtime_error("profilemover::deform(): curve network projection estimation failed.");
            return temp;
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Estimated cut-vert positions in " << elapsed.count() << " seconds"<< std::endl;

        // Compute per-face deformation matrix + assemble
        start = std::chrono::steady_clock::now();
        Eigen::MatrixXd y_h;
        int faceDef_success = CM.estimateFaceDeformations(y_h, CMheTohe, arap);
        if (faceDef_success != 1) {
            throw std::runtime_error("profilemover::deform(): face deformation estimate failed.");
            return temp;
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Computed face deformations in " << elapsed.count() << " seconds"<< std::endl;

        // Compute new positions
        start = std::chrono::steady_clock::now();
        Eigen::MatrixXd x_v = VtLV.solve(mVtL * (C * x_c - y_h));
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Solved for positions in " << elapsed.count() << " seconds"<< std::endl;
        
        start = std::chrono::steady_clock::now();
        std::vector<Eigen::Vector3d> final_pos = assembleFinalPositions(x_v, x_c);
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        // std::cout << "Assembled final positions in " << elapsed.count() << " seconds"<< std::endl;
        return final_pos;
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

    // Matrix forms

    int profilemover::assembleDiscreteCurvenetMats() {
        // TODO
        if ((f_dCN_flat.rows() != dCN.numHalfedges() && f_dCN_flat.cols() != 9) ||
            (x_dCN.rows() != dCN.numVerts() && x_dCN.cols() != 3) || 
            (f_dCN.rows() != 3*dCN.numHalfedges() && f_dCN.cols() != 3)) {
                f_dCN_flat.resize(dCN.numHalfedges(), 9);
                f_dCN_flat.setZero();
                x_dCN.resize(dCN.numVerts(), 3);
                x_dCN.setZero();
                f_dCN.resize(3*dCN.numHalfedges(), 3);
        }
        //dCN.computeDefGradOperators(f_dCN_flat, f_dCN, M_dCN_flat, M_3dCN_c);
        //dCN.computedCNVerts(x_dCN, M_dCN_c);
        return 1;
    }

    int profilemover::computeCDefGrads() {
        f_c = M_dCN_flat * f_dCN_flat;
        return 1;
    }

    int profilemover::computeFaceDefGrads() {
        Eigen::MatrixXd f_F_flat = M_v_F * f_v + M_c_F * f_c;
        return 1;
    }

    int profilemover::applyFaceDeformations() {
        // We don't need to fold the def grads together; we can just parallelize
        y_h = V.transpose() * x_v + C.transpose() * x_c;
        /*
        #pragma omp parallel for
        for (int he_idx = 0; he_idx < he_to_f.size(); he_idx++) {
            int he = he_to_f[he_idx].first;
            int f = he_to_f[he_idx].second;
            Eigen::Vector3d he_pos = y_h.row(he).transpose();
            // Use the transposed convention
            y_h(he, 0) = he_pos(0) * f_F_flat[f, 0] + he_pos(1) * f_F_flat[f, 3] + he_pos(2) * f_F_flat[f, 6];
            y_h(he, 1) = he_pos(0) * f_F_flat[f, 1] + he_pos(1) * f_F_flat[f, 4] + he_pos(2) * f_F_flat[f, 7];
            y_h(he, 2) = he_pos(0) * f_F_flat[f, 2] + he_pos(1) * f_F_flat[f, 5] + he_pos(2) * f_F_flat[f, 8];
        }
            */
        return 1;
    }

    int profilemover::computeCPositions() {
        x_c = M_dCN_c * x_dCN - M_3dCN_c * f_dCN * proj_c;
        return 1;
    }

    int profilemover::assembleFinalPositions(std::vector<Eigen::Vector3d>& newV) {
        Eigen::MatrixXd m = M_v_M * x_v + M_c_M * x_c;
        if (newV.size() != m.rows()) {
            newV.resize(m.rows());
            newV.clear();
        }
        #pragma omp parallel for
        for (int v = 0; v < m.rows(); v++) {
            newV[v] = m.row(v).transpose();
        }
        return 1;
    }

}   // namespace ProfileMover