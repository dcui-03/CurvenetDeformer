// mesh.hpp
#pragma once

#include "dcurvenet/dcurvenet.hpp"
#include "mesh_types.hpp"
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <vector>
#include <map>


namespace Mesh {

using vector2dList = std::vector<Eigen::Vector2d, Eigen::aligned_allocator<std::vector<Eigen::Vector2d>>>;

// Template mesh class that augments our standard mesh setup with some quantities we need
/* 
File Descriptors:
    - mesh_types.hpp: Structs for primal objects used by the mesh class (Verts, Edges, Faces, HalfEdges, etc.)
    - mesh.cpp: Mesh instance initialization functions and getters
    - mesh_utils.cpp: Standard mesh operations (projection, editing, etc.)
    - mesh_iter.cpp: Standard mesh iterators and queries (vertex-edge mapping, vertex umbrellas, edge-face mapping, etc.)
    - mesh_cut.cpp: Functions for cutting a mesh given a discrete curve network
*/
class mesh {
    public:
        // Constructor
        mesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List);

        // Add a discrete Curvenetwork Pointer
        bool applyDiscreteCurvenet(DCurvenet::dcurvenet* dCurvenet);
        // Assign a discrete curvenet index to a halfedge
        bool assignDCNtoHE(int he, const int dCN_idx, bool positive);
        // Add a reference mesh
        bool applyMeshRef(mesh* MRef);

        // Project a vertex onto the mesh
        // mesh_utils.cpp
        int computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, int& elIdx, bool snap = true, bool fast = true) const;
        
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

    private:
        // ------------- INITIALIZATION (mesh.cpp)  -----------------
        // Main init function
        bool initHalfEdgeMesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List);
        // Clear all mesh attributes
        bool clearMesh();
        // Internal function to precompute height functions on both planar/nonplanar faces
        Eigen::VectorXd computeFaceHeight(int f) const;
        // Internal function to precompute normals and areas on mesh structures
        double computeFVectorArea(int f, Eigen::Vector3d& fN);  // 1 face
        void computeFNormalsAreas();        // All faces
        // weight_fN weights by adjacent face areas
        double computeVNormalArea(int v, Eigen::Vector3d& vN, bool weight_fN = true);   // 1 vertex
        void computeVNormalsAreas(bool weight_fN = true);       // All vertices

        // Computes mean edge length on the mesh
        void computeMeanE();
        // Compute the length of the diagonal of the bounding box.
        void computeBBoxDiag();

        // ------------- UTILITIES (mesh_utils.cpp) -----------------

        // Inserts a vertex at a location into a data structure and appends to a corresponding face
        // NOTE: For safety, REQUIRE that vertex is attached to a real face
        // Returns the index of the new vertex
        int insertVertex(Eigen::Vector3d pos, int f, int dCN_idx = -1);
        // Topologically splits an existing edge by adding a new vertex.
        // NOTE: Added vertex does NOT need to lie on the edge
        // Returns index of the new vertex
        int splitEdge(int e, Eigen::Vector3d split_pos);
        // Main function that inserts a new edge connecting two vertices on a specified face
        // Returns the index of the new edge
        int insertEdge(int f, int v0, int v1, int dCN_idx0 = -1, int dCN_idx1 = -1, bool positive0 = true);
        // Helper that actually does the new edge insertion at an exact halfedge location
        int insertEdgeBetweenHEs(int f, int v0, int v1, int he0_prev, int v0_isolated, int he1_prev, int v1_isolated, int dCN_idx0, int dCN_idx1, bool positive0);
        // Helper that determines which halfedge to insert the new edge at.
        bool chooseEdgeInsertHE(int f, int v, int target, int& he_prev_out);

        // ------------- ITERATORS + QUERYING (mesh_iter.cpp) -----------------

        // Returns a CCW list of a vertex's OUTGOING halfedge indices
        std::vector<int> vertAdjHEs(int v) const;
        // Returns a CCW list of a vertex's adjacent vertices
        std::vector<int> vertAdjVerts(int v) const;
        // Returns a list of vertices in a loop from a given halfedge
        std::vector<int> mesh::vertLoop(int he) const;
        // Returns a CCW list of a vertex's adjacent faces
        std::vector<int> vertAdjFaces(int v) const;

        // Returns the endpoints of an edge in an arbitrary order.
        std::pair<int, int> edgeAdjVerts(int e) const;
        // Returns the adjacent face(s) of an edge (-1 indicates boundary)
        std::pair<int, int> edgeAdjFaces(int e) const;

        // Returns the halfedge index given the face index and edge index
        int halfedgeAtFaceEdge(int f, int e) const;
        // Get an arbitrary halfedge loop
        std::vector<int> halfedgeLoop(int he) const;

        // Returns a CCW list of a face's vertices
        std::vector<Eigen::Vector3d> faceAdjVerts(int f) const;
        std::vector<Eigen::Vector3d> mesh::faceAdjVerts(std::vector<int> fVerts) const;
        // Returns a CCW list of a face's vertex indices
        std::vector<int> faceAdjVertIdxs(int f) const;
        // Returns a CCW list of a face's half edges
        std::vector<int> faceAdjHalfEdges(int f) const;

        // Returns the outgoing boundary HE if a vertex is a boundary vertex, else returns -1
        int vertIsBoundary(int v, bool fast = true) const;
        // Returns true if halfedge is on boundary
        bool halfedgeIsBoundary(int he) const;

        // ------------- MESH CUTTING (mesh_cut.cpp) -----------------

        // Recursively computes straightest geodesic from a starting point given a starting direction, inserting new vertices as needed


        // ------------- ATTRIBUTES -----------------
        // List of primal mesh elements
        std::vector<Vert> V;
        std::vector<HalfEdge> HE;
        std::vector<Edge> E;
        std::vector<Face> F;
        // Counters for convenience
        int active_v = 0;
        int active_e = 0;
        int active_f = 0;
        // Vertex-to-Edge Map for easy indexing
        std::map<std::pair<int, int>, int> vertPairToHE;

        // Mean edge length on mesh
        double meanE;
        // AABB Diagonal length
        double bboxDiag;

        // Pointer to a dCN object if necessary
        DCurvenet::dcurvenet* dCN;
        bool dCN_initialized = false;

        // Pointer to a mesh object
        mesh* M_ref;
        bool M_ref_initialized = false;
};

}   // namespace Mesh