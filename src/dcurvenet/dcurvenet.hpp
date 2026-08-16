// dcurvenet.hpp
#pragma once

#include "dcurvenet_types.hpp"
#include "polynet/polynet.hpp"
#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
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

namespace Polynet {

// Discrete curvenet (i.e., polylines discretized from a curvenet)
class dcurvenet : public polynet {
    public:
        // NOTE: Vertices are copied directly from the Control of the curve network (CN),
        //       meaning they are in the same order and have the same corresponding indices.
        //       Curves are similar.
        // Takes the original curvenet and discretizes it
        // Alpha is the user-inputted sampling parameter
        dcurvenet(Curvenet::curvenet* CN, const Mesh::mesh* M = nullptr, bool sampleNaive = false);
        // Initialize with empty constructor
        dcurvenet();

        // --------- RUNTIME COMPUTATION -----------
        // Update with new curvenet positions and local frames
        void updateDiscCurveNet();

        // Propagate weights along curvenet
        int propagateWeights();

        // True iff the source curvenet has ordered (mesh-derived) connectivity at high-valence verts
        bool hasOrderedConnectivity() const;

        friend class Mesh::cutmesh;    // Friend class to access curvenet variables
        friend class ProfileMover::profilemover;
    protected:
        // No further inherited classes
    private:
        // --------- INITIALIZATION -----------
        // Add a vertex that matches an existing control
        int addVert(Curvenet::Control ctrl, int ctrl_idx = -1);
        // Add a vertex given its parameters (keeps vertData in lockstep with the inherited V)
        int addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n = Eigen::Vector3d::Zero());
        // Add a curve that matches an input curvenet curve
        int addCurve(int crv, bool sampleNaive = false);

        // --------- OTHER -----------
        // Check if a halfedge is the positive or negative side
        bool isPositiveHalfedge(int he) const;

        // Extra per-vertex/per-curve data linking back to the source curvenet, parallel to V/C
        std::vector<dCNVertData> vertData;
        std::vector<int> curveCNIdx;

        // Map input vertex index to local vertex index
        std::map<int, int> inputCtoV;
        // Map input curve index to local curve index
        std::map<int, int> inputCrvToC;

        // Pointer to parent curvenet
        Curvenet::curvenet* CN;
        // Which spline sampler was used to build this dcurvenet
        bool sampleNaive = false;
};

}   // namespace Polynet
