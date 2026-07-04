// cutmesh.hpp
#pragma once

#include "mesh/mesh.hpp"
#include "profilemover/profilemover.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "mesh/mesh_types.hpp"
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <vector>
#include <map>


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
        bool assignDCNtoHE(int he, const int dCN_idx, bool positive);
        // Add a reference mesh
        bool applyMeshRef(mesh* MRef);

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

        friend class ProfileMover::profilemover;
    protected:
    private:
        // ------------- INITIALIZATION (cutmesh.cpp)  -----------------
        bool copyFromMesh();

        // ------------- MESH CUTTING (cutmesh_cut.cpp) -----------------
        // Embed the pointed to discrete curve network into the mesh
        int embedCurves();
        int sortHalfEdges();
        int resetFaces();
        // Recursively computes straightest geodesic from a starting point and end point, inserting new vertices as needed
        // In order to check if the end is reached (i.e., when the start and end share a face), a "visibility" test is performed on the shared face via raycasting
        // Has a fast flag which simply checks for shared face(s) rather than doing visibility.
        int traceGeodesic(int start, int startF, int end, int endF, int dCN_idx0, int dCN_idx1, Eigen::Vector3d direc, bool reDirec, bool fast = true);
        // Actually cuts the mesh along the embedded vertices
        int cutMesh();

        // ------------- UTILITIES (mesh_utils.cpp) -----------------

        // Inserts a vertex at a location into a data structure
        // Returns the index of the new vertex
        // TODO: Move to protected part of mesh class
        int insertVertex(Eigen::Vector3d pos,
                         Eigen::Vector3d n, 
                         int label = 0, 
                         int dCN_idx = -1, 
                         int ref_Type = 0, 
                         int ref_Idx = -1, 
                         Eigen::Vector3d proj = Eigen::Vector3d::Zero());
        int insertVertex(Vert splitV);
        Vert createVertex(Eigen::Vector3d pos,
                          Eigen::Vector3d n, 
                          int label = 0, 
                          int dCN_idx = -1, 
                          int ref_Type = 0, 
                          int ref_Idx = -1, 
                          Eigen::Vector3d proj = Eigen::Vector3d::Zero());
        // Topologically split an edge with an existing vertex
        int splitEdge(int e, int new_v);
        // Insert an edge between two existing vertices
        int insertEdge(int v0, int v1, int dCN_idx0, int dCN_idx1);

        // Pointer to a reference Mesh object
        mesh* M;
        // Pointer to a dCN object if necessary
        DCurvenet::dcurvenet* dCN;
        // Map from dCN vert indices to cut-mesh vertex indices
        std::map<int, std::vector<int>> dCNVtoV;
        // Map from dCN halfedge indices to cut-mesh halfedge indices
        std::map<int, std::vector<int>> dCNHEtoHE;
};

}   // namespace Mesh