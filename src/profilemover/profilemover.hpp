// profilemover.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "mesh/mesh.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>


namespace ProfileMover {

class profilemover {
    public:
        // Constructor, which first builds the mesh
        profilemover(std::vector<Eigen::Vector3d>& meshV, std::vector<std::vector<int>>& meshF);
        profilemover();

        // Precompute cut-mesh and operators
        // Takes as input the necessary items to construct the curve network
        void precomputation(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents, std::vector<std::array<int, 4>> Splines, int alpha = 5);

        // Apply deformation given the new control and tangent locations (connectivity should be same)
        // Returns new mesh positions as an Nx3 matrix
        void deform(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents);
    protected:
        // No class inheritance
    private:
        // Computes deformation matrix, where each row is the flattened deformation gradient of a dCN halfedge
        Eigen::MatrixXd computeFlatDefGrads();

        // TODO: Intermediary stages
        // Per-cutmesh per-face deformation gradient
        // Deformed projection needed for second opt
        Eigen::MatrixXd estimateProjectionDefs();

        // Store copy of cut mesh
        Mesh::mesh M;

        // NOTE: no need to store updated states curvenet and dcurvenet, we have to create new copies at execution time
        // Store the neutral curvenet (spline)
        // TODO: Do we need the original curvenet? --> Only for resets
        // We do need the original nDCurvenet so that we can compute deformation gradients
        Curvenet::curvenet CN;
        // Store the neutral discrete curvenet
        DCurvenet::dcurvenet dCN;
        bool dCN_init = false;
        // Store neutral cut-mesh
        //CutMesh::cutmesh nCutmesh;

        
        // Store operators
        // TODO: Need functions to compute V and C
        // TODO: Instead of storing VtLV, store its factorization.
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> VtLV;
        Eigen::SparseMatrix<double> mVtL;
        Eigen::SparseMatrix<double> V;
        Eigen::SparseMatrix<double> C;
        
};

}   // namespace ProfileMover