// polynet.hpp
#pragma once

#include "polynet_types.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>

namespace Polynet {

// Generic polyline network (topology, mesh binding, weight diffusion)
/*
File Descriptors:
    - polynet_types.hpp: Structs for primal objects used by the polynet class (Verts, Edges, Curves, HalfEdges, etc.)
    - polynet.cpp: Construction, topology helpers, weight propagation
    - polynet_trace.cpp: Trace out curves at init
*/
class polynet {
    public:
        // Builds its own curves from a raw vertex/edge soup
        polynet(const std::vector<Eigen::Vector3d>& V_list, const std::vector<std::array<int, 2>>& E_list, const Mesh::mesh* M = nullptr);
        // Initialize with empty constructor
        polynet();

        // --------- GETTERS -----------
        const int numVerts() const { return V.size(); }
        const int numHalfedges() const { return HE.size(); }

        // --------- RUNTIME COMPUTATION -----------
        // Propagate weights along the network using V[].fixed_w/.w as the boundary condition
        int propagateWeights();

    protected:
        // --------- INITIALIZATION -----------
        // Add a vertex given its parameters
        int addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n = Eigen::Vector3d::Zero());
        // Add an edge and return the index of the new edge
        int addEdge(int origin, int dest, int prev_he0 = -1, int next_he1 = -1, int c = -1);
        // Rewire incoming/outgoing halfedges of a vertex such that topology is correct
        // At valence 1/2 this is always unambiguous. At valence 3+, prev/next are only
        // wired if `ordered` is true (i.e. adjHE is a genuine CCW order); otherwise they
        // are left at -1, marking the halfedge chain as ending there.
        int rewireVertAdjHE(int v, bool ordered);
        // Compute projection data
        int computeProjData(const Mesh::mesh* M);

        // --------- VALENCE-BASED CURVE TRACING (polynet_trace.cpp) -----------
        // Compute what kind of vertex each vert is using the valence of adjacent edges
        int assignVertType(int v);
        int assignVertTypeAll();
        // Trace out curves
        int traceCurves();
        // Helper to find next traced edge
        int nextHEFromVert(int curr_he, int curr_end);

        // Attributes as lists
        std::vector<Vert> V;
        std::vector<HalfEdge> HE;
        std::vector<Edge> E;
        std::vector<Curve> C;

        // True once a mesh has been bound (at construction or later)
        bool setMesh = false;
};

}   // namespace Polynet
