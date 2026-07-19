// cutmesh.hpp
#pragma once

#include "mesh/mesh.hpp"
#include "mesh/mesh_types.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

namespace DCurvenet {
    class dcurvenet;
}

namespace Mesh {

// Cutmesh class that augments our standard mesh setup with some quantities we need
/* 
File Descriptors:
    - mesh_types.hpp: Structs for primal objects used by the mesh class (Verts, Edges, Faces, HalfEdges, etc.)
    - mesh.cpp: Mesh instance initialization functions and getters
    - mesh_utils.cpp: Standard mesh operations (projection, editing, etc.)
    - mesh_iter.cpp: Standard mesh iterators and queries (vertex-edge mapping, vertex umbrellas, edge-face mapping, etc.)
    - mesh_cut.cpp: Functions for cutting a mesh given a discrete curve network
*/
class cutmesh : public mesh {
    public:
        // Constructor
        // Copies in vertex and edge data from the reference mesh then applies the dCN to embed the curves
        cutmesh(mesh* MRef, DCurvenet::dcurvenet* dCN);
        cutmesh();

        // Add a discrete Curvenetwork Pointer
        bool applyDiscreteCurvenet(DCurvenet::dcurvenet* dCurvenet);
        // Assign a discrete curvenet index to a halfedge
        bool assignDCNtoHE(int he, const int dCN_idx);
        // Add a reference mesh
        bool applyMeshRef(mesh* MRef);

        // Helpers for profilemover class
        void computeHEMap(std::vector<int>& heToCMhe, std::map<int, int>& CMheTohe);
        int computeVMatrix(Eigen::SparseMatrix<double>& V_mat, std::vector<int>& vToCM, std::map<int, int>& mToV, const std::vector<int>& heToCMhe);
        int computeCMatrix(Eigen::SparseMatrix<double>& C_mat, std::vector<int>& cToCM, std::map<int, std::vector<int>>& mToC, const std::vector<int>& heToCMhe);
        
        int compute_dCNMaps(Eigen::SparseMatrix<double>& M_dCN_flat,
                             Eigen::SparseMatrix<double>& M_dCN_c,
                             const std::vector<int>& cToCM);
        int computeProjMatrix(Eigen::MatrixXd& projVecs, const std::vector<int>& cToCM);
        int computeFaceOps(Eigen::SparseMatrix<double>& M_v_F,
                            Eigen::SparseMatrix<double>& M_c_F,
                            std::vector<int>& M_he_F,
                            Eigen::MatrixXd& x_h,
                            const std::vector<int>& vToCM,
                            const std::vector<int>& cToCM,
                            const std::vector<int>& heToCMhe,
                            const std::map<int, int>& CMheTohe);
        int computeAssemblyOps(Eigen::SparseMatrix<double>& M_v_M,
                               Eigen::SparseMatrix<double>& M_c_M,
                               const std::map<int, int>& mToV,
                               const std::map<int, std::vector<int>>& mToC,
                               int num_M, int num_V, int num_C);
        
        int computeHELaplacian(Eigen::SparseMatrix<double>& L, std::map<int, int>& CMheTohe);
        // Compute the deformation gradient on the initial cutmesh dCN verts
        int computeDefGrads(Eigen::MatrixXd& defGrads, const std::vector<int>& cToCM);
        // Apply solved deformation gradients to the cutmesh
        void applyDefGrads(const Eigen::MatrixXd& defGrads, const std::vector<int>& vToCM);
        // Estimate projected curvenet positions
        int estimateCNPositions(Eigen::MatrixXd& cnPos, const std::vector<int>& cToCM);
        // Estimate the deformed faces
        int estimateFaceDeformations(Eigen::MatrixXd& faceDef, const std::map<int, int>& CMheTohe, bool arap = false);

        // POLYSCOPE reformatting
        int polyscopeFormat(Eigen::MatrixXd& Verts, std::vector<std::vector<int>>& Faces, 
                            std::vector<glm::vec3>& VertN, std::vector<glm::vec3>& FaceN, 
                            std::vector<glm::vec3>& cornerIdx,
                            std::vector<glm::vec3>& projVecs) const;

        // Getters
        Eigen::Vector3d getVPos(int v) const;
        Eigen::Vector3d getNormal(int elType, int elIdx) const;
        Eigen::Vector3d getVNormal(int v) const;
        Eigen::Vector3d getENormal(int e) const;
        Eigen::Vector3d getFNormal(int f) const;

        // Get mean edge length
        double getMeanE() const;
        // Get bbox diagonal length
        double getBBoxDiag() const;

        // friend class ProfileMover::profilemover;
    protected:
    private:
        // ------------- INITIALIZATION (cutmesh.cpp)  -----------------
        bool copyFromMesh();

        // ------------- MESH CUTTING (cutmesh_cut.cpp) -----------------
        // Embed the pointed to discrete curve network into the mesh
        int embedCurves();
        // Sort all halfedges so that their next/prev are correct
        int sortHalfEdges();
        // Reset faces
        int resetFaces();
        // Find and deactivate any projected curves/loops that are not attached to the original mesh
        int deactivateIsolatedCuts();
        // Actually cuts the mesh along the embedded vertices
        int cutMesh();

        // ------------- UTILITIES (mesh_utils.cpp) -----------------

        // Inserts a vertex at a location into a data structure
        // Returns the index of the new vertex
        // TODO: Move to protected part of mesh class
        int insertVertex(Eigen::Vector3d pos,
                         Eigen::Vector3d n, 
                         int label = 0, 
                         int cornerIdx = -1, 
                         int ref_Type = 0, 
                         int ref_Idx = -1, 
                         Eigen::Vector3d proj = Eigen::Vector3d::Zero(),
                         Eigen::Matrix3d defGrad = Eigen::Matrix3d::Identity());
        int insertVertex(Eigen::Vector3d pos, 
                        Eigen::Vector3d n, 
                        int label = 0, 
                        int cornerIdx = -1, 
                        vertProjData projData = vertProjData({-1, -1}), 
                        vertDeformData defData = vertDeformData({Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity()}));
        int insertVertex(Vert splitV);
        // Topologically split an edge with an existing vertex
        int splitEdge(int e, int new_v);
        // Insert an edge between two existing vertices
        int insertEdge(int v0, int v1, int dCN_idx0, int dCN_idx1);

        // Pointer to a reference Mesh object
        mesh* M;
        // Pointer to a dCN object if necessary
        DCurvenet::dcurvenet* dCN;
};

}   // namespace Mesh