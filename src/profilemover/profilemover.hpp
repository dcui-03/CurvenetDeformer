// profilemover.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "mesh/mesh.hpp"
#include "cutmesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include <vector>
#include <array>
#include <map>
#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>


namespace ProfileMover {

class profilemover {
    public:
        // Constructor, which first builds the mesh
        profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF, 
                               const std::vector<Eigen::Vector3d> Controls, const std::vector<Eigen::Vector3d> Tangents, 
                               const std::vector<std::array<int, 4>> Splines, int alpha = 5, bool arap = false);
        // Only apply mesh
        profilemover(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF);
        profilemover();

        // Getters in case we need it
        const Mesh::mesh& mesh() const;
        const Mesh::cutmesh& cutmesh() const;
        const Curvenet::curvenet& curvenet() const;
        const DCurvenet::dcurvenet& discreteCurvenet() const;
        
        void toggleARAP(bool toggle);
        void applyMesh(const std::vector<Eigen::Vector3d>& meshV, const std::vector<std::vector<int>>& meshF);
        void applyCurvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, int alpha = 5);
        void computeDiscreteCurvenet();
        void computeCutMesh();
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

        // Matrix-forms of runtime computation
        int assembleDiscreteCurvenetMats();
        int computeCDefGrads();
        int computeFaceDefGrads();
        int applyFaceDeformations();
        int computeCPositions();
        int assembleFinalPositions(std::vector<Eigen::Vector3d>& newV);

        // ARAP-style deformations
        bool arap = false;
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
        // Map from halfedge indices in V and C to the cutmesh
        std::vector<int> heToCMhe;
        std::map<int, int> CMheTohe;

        // Map from the orginal mesh's vertices to the C and V vertices
        // Note that I'm opting to use a map here, since we need to split between C and V
        // Instead of maps, store these as precomputed operators?
        std::map<int, std::vector<int>> mToC;
        std::map<int, int> mToV;

        // Matrix forms of the relevant cutmesh components
        Eigen::MatrixXd proj_c;                         // Dense matrix with the proj vectors for each C cut-vertex
        Eigen::SparseMatrix<double> he_pos;             // Sparse, block matrix of corner positions per face
        Eigen::MatrixXd f_dCN_flat;                     // Dense matrix of flattened dCN-based deformation gradients
        Eigen::MatrixXd f_dCN;                          // Dense matrix of compressed dCN-based deformation gradients
        Eigen::MatrixXd f_faces;                        // Dense matrix of def grads on cut-faces
        Eigen::MatrixXd f_c;                            // Dense matrix of cutmesh C cut-vertex flattened def grads
        Eigen::MatrixXd f_v;                            // Dense matrix of cutmesh V cut-vertex flattened def grads
        Eigen::MatrixXd x_c;                            // Dense matrix of cutmesh C cut-vertex estimated positions
        Eigen::MatrixXd x_v;                            // Dense matrix of cutmesh V cut-vertex solved positions
        Eigen::MatrixXd x_dCN;                          // Dense matrix of dCN vertices
        Eigen::MatrixXd y_h;

        // Sparse maps
        Eigen::SparseMatrix<double> M_dCN_flat;         // Sparse averaging map from dCN halfedges to cutmesh C
        Eigen::SparseMatrix<double> M_v_F;              // Sparse matrix from V cut-vertices to faces
        Eigen::SparseMatrix<double> M_c_F;              // Sparse matrix from C cut-vertices to faces
        Eigen::SparseMatrix<double> M_v_M;              // Sparse matrix from V cut-vertices to the mesh
        Eigen::SparseMatrix<double> M_c_M;              // Sparse matrix from C cut-vertices to the mesh
        Eigen::SparseMatrix<double> M_dCN_c;            // Sparse matrix from dCN verts to cutmesh C
        Eigen::SparseMatrix<double> M_3dCN_c;           // Sparse matrix from dCN def grads to cutmesh C
        std::vector<std::pair<int, int>> he_to_f;       // Map from local halfedge indices to their corresponding local face indices

        
        // Store operators
        // TODO: Need functions to compute V and C
        // TODO: Instead of storing VtLV, store its factorization.
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> VtLV;
        Eigen::SparseMatrix<double> mVtL;
        Eigen::SparseMatrix<double> V;
        Eigen::SparseMatrix<double> C;
        
};

}   // namespace ProfileMover