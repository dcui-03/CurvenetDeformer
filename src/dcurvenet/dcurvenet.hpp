// dcurvenet.hpp
#pragma once

#include "dcurvenet_types.hpp"
#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <map>

namespace Mesh {
    class cutmesh;
}

namespace Curvenet {
    class curvenet;
}

namespace ProfileMover {
    class profilemover;
}

namespace DCurvenet {

// Discrete curvenet (i.e., polylines)
class dcurvenet {
    public:
        // NOTE: Vertices are copied directly from the Control of the curve network (CN),
        //       meaning they are in the same order and have the same corresponding indices.
        //       Curves are similar.
        // Takes the original curvenet and discretizes it
        // Alpha is the user-inputted sampling parameter
        dcurvenet(Curvenet::curvenet* CN, Mesh::mesh* M);
        // Initialize with empty constructor
        dcurvenet();

        // --------- GETTERS -----------
        const int numVerts() const { return V.size(); }
        const int numHalfedges() const { return HE.size(); }

        // --------- RUNTIME COMPUTATION -----------
        // Update with new curvenet positions and local frames
        void updateDiscCurveNet();

        // Propagate weights along curvenet
        int propagateWeights();

        friend class Mesh::cutmesh;    // Friend class to access curvenet variables
        friend class ProfileMover::profilemover;
    protected:
        // No inherited classes
    private:
        // --------- INITIALIZATION -----------
        // Add a vertex
        int addVert(Curvenet::Control ctrl, int ctrl_idx = -1);
        int addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n = Eigen::Vector3d::Zero(), int ctrl_idx = -1, int ctrl_type = -1, int adjSize = 0);
        // Add an edge and return the index of the new edge
        int addEdge(int origin, int dest, int prev_he0 = -1, int next_he1 = -1, int c = -1);
        // Add a curve
        int addCurve(int crv);
        // Rewire incoming/outgoing halfedges of an intersection vertex such that topology is correct
        int rewireVertAdjHE(int v);
        // Compute projection data
        int computeProjData(Mesh::mesh* M);

        // --------- OTHER -----------
        // Check if a halfedge is the positive or negative side
        bool isPositiveHalfedge(int he) const;

        // Attributes as lists
        std::vector<Vert> V;
        std::vector<HalfEdge> HE;
        std::vector<Edge> E;
        std::vector<Curve> C;
        
        // Map input vertex index to local vertex index
        std::map<int, int> inputCtoV;
        // Map input curve index to local curve index
        std::map<int, int> inputCrvToC;

        // Pointer to parent curvenet
        Curvenet::curvenet* CN;
};

}   // namespace DCurvenet