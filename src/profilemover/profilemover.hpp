// profilemover.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "mesh/mesh.hpp"
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>


namespace ProfileMover {

class profilemover {
    public:
        // Constructor, which first builds the mesh
        profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF, 
                               const std::vector<Eigen::Vector3d> Controls, const std::vector<Eigen::Vector3d> Tangents, 
                               const std::vector<std::array<int, 4>> Splines, int alpha = 5);
        profilemover();
        
        void applyMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF);
        void applyCurvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, int alpha = 5);

        // Apply deformation given the new control and tangent locations (connectivity should be same)
        // Returns new mesh positions as an Nx3 matrix
        std::vector<Eigen::Vector3d> deform(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents);
    protected:
        // No class inheritance
    private:
        // Precompute cut-mesh and operators
        // Takes as input the necessary items to construct the curve network
        void precomputation();
        // Assemble final positions into our standard data type
        std::vector<Eigen::Vector3d> assembleFinalPositions(Eigen::MatrixXd x_v, Eigen::MatrixXd x_c);

        // Store copy of mesh
        Mesh::mesh M;
        bool M_init = false;
        // Store the curvenet
        Curvenet::curvenet CN;
        bool CN_init = false;
        // Store the neutral discrete curvenet
        DCurvenet::dcurvenet dCN;
        bool dCN_init = false;
        // Store the cut-mesh
        Mesh::cutmesh CM;
        bool CM_init = false;

        // Map from matrices (C, V) to cutmesh cut-vertices and vice versa
        // BE CAREFUL about inactive vertices
        // NOTE: Could use a std::map insted, but this is more intuitive, since we will end up parallelizing these steps
        std::vector<int> vToCM;
        // Map from matrices to cutmesh regular vertices and vice versa
        std::vector<int> cToCM;

        // Map C to the associated mesh vertices if they land on vertices
        std::vector<std::vector<int>> mToC;
        std::vector<int> mToV;
        
        // Store operators
        // TODO: Need functions to compute V and C
        // TODO: Instead of storing VtLV, store its factorization.
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> VtLV;
        Eigen::SparseMatrix<double> mVtL;
        Eigen::SparseMatrix<double> V;
        Eigen::SparseMatrix<double> C;
        
};

}   // namespace ProfileMover