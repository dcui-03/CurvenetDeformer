#include "profilemover.hpp"

#include "mesh/mesh.hpp"
#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
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
        
        //precomputation();
        precomputeOps();
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

    // Weights
    void profilemover::assignWeight(int cnVert, bool fixed, double w) {
        if (!CN_init || !dCN_init) {
            return;
        }
        int success = CN.assignWeight(cnVert, fixed, w);
        if (success != 1) {
            return;
        } 
        // Rebuild weight operator
        success = dCN.propagateWeights();
        if (success != 1) {
            return;
        }
        recompute_weights = true;
        return;
    }
    void profilemover::clearWeights() {
        if (!CN_init || !dCN_init) {
            return;
        }
        CN.resetWeights();
        // Rebuild weight operator
        int success = dCN.propagateWeights();
        if (success != 1) {
            return;
        }
        recompute_weights = true;
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
            dCN = DCurvenet::dcurvenet(&CN, &M);
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
        int C_success = CM.computeCMatrix(C, c_rest, cToCM, mToC, heToCMhe);
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
    std::vector<Eigen::Vector3d> profilemover::deform(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents) {
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
    void profilemover::precomputeOps() {
        if (!M_init || !CN_init || !dCN_init || !CM_init) {
            throw std::runtime_error("profilemover::precomputeOps(): curvenet, mesh, cutmesh, or discrete curvenet not initialized.");
            return;
        }
        auto Op_start = std::chrono::steady_clock::now();
        // Compute operators
        CM.computeHEMap(heToCMhe, CMheTohe);
        int V_success = CM.computeVMatrix(V, vToCM, mToV, heToCMhe);
        std::cout << "Computed V" << std::endl;
        if (V_success != 1) {
            throw std::runtime_error("profilemover::precomputateOps(): Unable to construct matrix V.");
            return;
        }
        int C_success = CM.computeCMatrix(C, c_rest, cToCM, mToC, heToCMhe);
        std::cout << "Computed C" << std::endl;
        if (C_success != 1) {
            throw std::runtime_error("profilemover::precomputateOps(): Unable to construct matrix C.");
            return;
        }
        // Pre-compute a matrix of projection vectors
        int proj_success = CM.computeProjMatrix(proj_c, cToCM);
        std::cout << "Computed proj Matrix" << std::endl;
        if (proj_success != 1) {
            throw std::runtime_error("profilemover::precomputateOps(): Unable to compute projection matrix.");
            return;
        }

        int dcn_success = CM.compute_dCNMaps(M_dCN_flat, M_dCN_c, cToCM);
        std::cout << "Computed dCN maps" << std::endl;
        if (dcn_success != 1) {
            throw std::runtime_error("profilemover::precomputeOps(): Unable to compute dCN operators.");
        }
        int weights_success = CM.computeWeightOps(weights, cToCM);
        std::cout << "Computed weight operator" << std::endl;
        if (weights_success != 1) {
            throw std::runtime_error("profilemover::precomputeOps(): Unable to compute weight operator.");
        }
        recompute_weights = false;

        int face_success = CM.computeFaceOps(M_v_F, M_c_F, M_he_F, x_h, vToCM, cToCM, heToCMhe, CMheTohe);
        std::cout << "Computed face ops" << std::endl;
        if (face_success != 1) {
            throw std::runtime_error("profilemover::precomputeOps(): Unable to compute face operators.");
        }

        int assembly_success = CM.computeAssemblyOps(M_v_M, M_c_M, mToV, mToC, M.getNumActiveV(), vToCM.size(), cToCM.size());
        std::cout << "Computed assembly ops" << std::endl;
        if (assembly_success != 1) {
            throw std::runtime_error("profilemover::precomputeOps(): Unable to compute final assembly operators.");
        }
        // Face-based Laplacian with values pushed onto halfedges
        Eigen::SparseMatrix<double> L;
        int L_success = CM.computeHELaplacian(L, CMheTohe);
        std::cout << "Computed laplacian" << std::endl;
        if (L_success != 1) {
            throw std::runtime_error("profilemover::precomputateOps(): Error constructing halfedge Laplacian.");
            return;
        }
        // mVtL
        mVtL = -1 * V.transpose() * L;
        // Factor V^TLV
        Eigen::SparseMatrix<double> VtLV_Mat = V.transpose() * L * V;
        // VtLV_Mat = V.transpose() * L * V;
        VtLV.analyzePattern(VtLV_Mat);
        VtLV.factorize(VtLV_Mat);
        // Iterative solver with preconditioner
        // VtLV.compute(VtLV_Mat);
        // VtLV.setMaxIterations(20);
        // VtLV.setTolerance(0.001);
        auto Op_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> Op_elapsed = Op_end - Op_start;
        std::cout << "Operators computed in time " << Op_elapsed.count() << " seconds" << std::endl;

        f_v.resize(V.cols(), 9);
        f_v.setZero();

        f_c.resize(C.cols(), 9);
        f_c.setZero();

        x_v.resize(V.cols(), 3);
        x_v.setZero();

        x_c.resize(C.cols(), 3);
        x_c.setZero();
        return;
    }

    std::vector<Eigen::Vector3d> profilemover::deformOps(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents) {
        std::vector<Eigen::Vector3d> newV;
        if (!M_init || !CN_init || !dCN_init || !CM_init) {
            throw std::runtime_error("profilemover::deformOps(): mesh, curvenet, discrete curvenet, or cutmesh not computed");
            return newV;
        }
        // 1. Compute new curvenet
        auto start = std::chrono::steady_clock::now();
        CN.updateCurveNet(Controls, Tangents);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Updated Curvenet in " << elapsed.count() << " seconds"<< std::endl;
        // 2. Compute new discrete curvenet and frames
        start = std::chrono::steady_clock::now();
        dCN.updateDiscCurveNet();
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        std::cout << "Updated discrete Curvenet in " << elapsed.count() << " seconds"<< std::endl;
        if (recompute_weights) {
            CM.computeWeightOps(weights, cToCM);
            recompute_weights = false;
        }
        // FIRST SOLVE: Deformation gradients
        // Compute flattened deformation gradient matrix
        start = std::chrono::steady_clock::now();
        int DG_success = assembleDiscreteCurvenetMats();
        if (DG_success != 1) {
            throw std::runtime_error("profilemover::deformOps(): deformation gradient computation failed.");
            return newV;
        }
        // Constraint deformation gradients
        f_c = M_dCN_flat * f_dCN_flat;
        if (weight_defgrads) {
            // Multiply by weights
            weightDefGrads();
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        std::cout << "Computed cut-vert def grads in " << elapsed.count() << " seconds"<< std::endl;
        // Solve system to get interpolated def grads
        start = std::chrono::steady_clock::now();
        f_v = VtLV.solve(mVtL * (C * f_c));
        // Eigen::MatrixXd rhs_f = mVtL * (C * f_c);
        // f_v = VtLV.solveWithGuess(rhs_f, f_v);
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        std::cout << "Computed all def grads in " << elapsed.count() << " seconds"<< std::endl;

        // SECOND SOLVE: Positions
        // Estimate new projected positions using the distributed def grads
        start = std::chrono::steady_clock::now();
        int cnPos_success = applyDefGradsToProj();
        if (cnPos_success != 1) {
            throw std::runtime_error("profilemover::deformOps(): curve network projection estimation failed.");
            return newV;
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        std::cout << "Applied def grads to projected vectors in " << elapsed.count() << " seconds"<< std::endl;

        start = std::chrono::steady_clock::now();
        if (weight_pos) {
            Eigen::MatrixXd q_new  = M_dCN_c * x_dCN;
            Eigen::MatrixXd q_rest = M_dCN_c * c_rest;
            for (int c = 0; c < q_new.rows(); c++) {
                double w = std::clamp(weights[c], 0.0, 1.0);
                q_new.row(c) = q_rest.row(c) + w * (q_new.row(c) - q_rest.row(c));
            }
            x_c = q_new - x_c;
        } else {
            x_c = M_dCN_c * x_dCN - x_c;
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        std::cout << "Estimated cut-vert positions in " << elapsed.count() << " seconds"<< std::endl;

        // Compute per-face deformation matrix + assemble
        start = std::chrono::steady_clock::now();
        Eigen::MatrixXd f_F_flat = M_v_F * f_v + M_c_F * f_c;
        int faceDef_success = applyFaceDeformations(f_F_flat);
        if (faceDef_success != 1) {
            throw std::runtime_error("profilemover::deformOps(): face deformation estimate failed.");
            return newV;
        }
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
        std::cout << "Computed face deformations in " << elapsed.count() << " seconds"<< std::endl;

        // Compute new positions
        start = std::chrono::steady_clock::now();
        // x_v = VtLV.solveWithGuess(mVtL * (C * x_c - y_h), x_v);
        x_v = VtLV.solve(mVtL * (C * x_c - y_h));
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
         std::cout << "Solved for positions in " << elapsed.count() << " seconds"<< std::endl;
        
        start = std::chrono::steady_clock::now();
        int final_success = assembleFinalPositions(newV);
        end = std::chrono::steady_clock::now();
        elapsed = end - start;
         std::cout << "Assembled final positions in " << elapsed.count() << " seconds"<< std::endl;
        return newV;
    }

    int profilemover::assembleDiscreteCurvenetMats() {
        if (f_dCN_flat.rows() != dCN.numHalfedges() || f_dCN_flat.cols() != 9) {
            f_dCN_flat.resize(dCN.numHalfedges(), 9);
        }

        if (x_dCN.rows() != dCN.numVerts() || x_dCN.cols() != 3) {
            x_dCN.resize(dCN.numVerts(), 3);
        }

        return dCN.computedCNMats(f_dCN_flat, x_dCN);
    }

    int profilemover::weightDefGrads() {
        if (weights.rows() != f_c.rows()) {
            return -1;
        }
        Eigen::VectorXd Id_flat = Utils::flattenMatrix3d(Eigen::Matrix3d::Identity());
        #pragma omp parallel for
        for (int r = 0; r < f_c.rows(); r++) {
            f_c.row(r) = (1 - weights[r]) * Id_flat.transpose() + weights[r] * f_c.row(r);
        }
        return 1;
    }

    int profilemover::applyDefGradsToProj() {
        int num_C = cToCM.size();
        x_c.resize(num_C, 3);

        #pragma omp parallel for
        for (int c = 0; c < num_C; c++) {
            Eigen::VectorXd f = f_c.row(c).transpose();
            double px = proj_c(c, 0);
            double py = proj_c(c, 1);
            double pz = proj_c(c, 2);

            x_c(c, 0) = f_c(c, 0) * px + f_c(c, 3) * py + f_c(c, 6) * pz;
            x_c(c, 1) = f_c(c, 1) * px + f_c(c, 4) * py + f_c(c, 7) * pz;
            x_c(c, 2) = f_c(c, 2) * px + f_c(c, 5) * py + f_c(c, 8) * pz;
        }
        return 1;
    }

    int profilemover::applyFaceDeformations(const Eigen::MatrixXd& f_F_flat) {
        const int num_h = static_cast<int>(heToCMhe.size());
        y_h.resize(num_h, 3);
        // Non-ARAP version
        if (!arap) {
            #pragma omp parallel for
            for (int he = 0; he < num_h; he++) {
                int f = M_he_F[he];

                double x = x_h(he, 0);
                double y = x_h(he, 1);
                double z = x_h(he, 2);
                // Compute the deformation directly (no folding)
                y_h(he, 0) = x * f_F_flat(f, 0) + y * f_F_flat(f, 3) + z * f_F_flat(f, 6);
                y_h(he, 1) = x * f_F_flat(f, 1) + y * f_F_flat(f, 4) + z * f_F_flat(f, 7);
                y_h(he, 2) = x * f_F_flat(f, 2) + y * f_F_flat(f, 5) + z * f_F_flat(f, 8);
            }
            return 1;
        }

        // ARAP branch only does polar decomposition
        std::vector<Eigen::Matrix3d> faceTransform(f_F_flat.rows());

        #pragma omp parallel for
        for (int f = 0; f < f_F_flat.rows(); f++) {
            Eigen::Matrix3d F = Utils::compressVector9d(f_F_flat.row(f).transpose());
            Eigen::Matrix3d R, S;
            Utils::polarDecomposition(F, R, S);
            faceTransform[f] = R;
        }

        #pragma omp parallel for
        for (int he = 0; he < num_h; he++) {
            int f = M_he_F[he];
            y_h.row(he) = x_h.row(he) * faceTransform[f].transpose();
        }

        return 1;
    }

    int profilemover::assembleFinalPositions(std::vector<Eigen::Vector3d>& newV) {
        Eigen::MatrixXd m = M_v_M * x_v + M_c_M * x_c;
        if (newV.size() != m.rows()) {
            newV.resize(m.rows());
        }
        #pragma omp parallel for
        for (int v = 0; v < m.rows(); v++) {
            newV[v] = m.row(v).transpose();
        }
        return 1;
    }

}   // namespace ProfileMover