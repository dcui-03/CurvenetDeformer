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
        profilemover(std::vector<Eigen::Vector3d>& meshV, std::vector<std::vector<int>>& meshF, 
                    std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents, std::vector<std::array<int, 4>> Splines, int alpha = 5);
        profilemover();

        // Precompute cut-mesh and operators
        // Takes as input the necessary items to construct the curve network
        void precomputation();

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

        // Store copy of mesh
        Mesh::mesh M;
        // Store a copy of cut-mesh
        Mesh::cutmesh CM;
        
        Curvenet::curvenet CN;
        // Store the neutral discrete curvenet
        DCurvenet::dcurvenet dCN;
        bool dCN_init = false;

        // Map from matrices (C, V) to cutmesh cut-vertices and vice versa
        // BE CAREFUL about inactive vertices
        // NOTE: Could use a std::map insted, but this is more intuitive, since we will end up parallelizing these steps
        std::vector<std::pair<int, int>> cutmeshVtoC_idx;
        // Map from matrices to cutmesh regular vertices and vice versa
        std::vector<std::pair<int, int>> cutmeshVtoV_idx;

        // Note: We also may need some mapping from cut-verts w/ label 1 and coincident to a mesh vertex in order to assemble final positions
        //       BUT this may be best assembled in the cutmesh class rather than in the profilemover class.
        
        // Store operators
        // TODO: Need functions to compute V and C
        // TODO: Instead of storing VtLV, store its factorization.
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> VtLV;
        Eigen::SparseMatrix<double> mVtL;
        Eigen::SparseMatrix<double> V;
        Eigen::SparseMatrix<double> C;
        
};

}   // namespace ProfileMover