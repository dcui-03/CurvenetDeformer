// curvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

// File with basic structs used by mesh class

namespace Curvenet {
    // Control vertices
    // NOTE: Each control stores an (ordered) list of outgoing halfedges
    //       This is in case a spline starts and ends at the same control
    struct Control {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        std::vector<int> adjHE;        // Outgoing Halfedge list
        bool active = true;     // For safety, say if the component is active (ignore for now)
        bool sorted = false;    // Safety flag. True when outgoing halfedges are sorted
    };

    // NOTE: Halfedge iteration indices (next, prev) prioritize easy iteration over their associated controls
    //       As a result, they may not be ideal for standard halfedge traversal tasks.
    struct HalfEdge {
        int origin = -1;  // To control
        int s = -1;     // associated spline
        // Iterators: Note that if we are at an endpoint, leave entry as -1
        int twin = -1;
        int next = -1;
        int prev = -1;
        Eigen::Vector3d tan;   // Tangent vector (defined in global coordinates, NOT relative to control)

        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Splines for easy iteration
    struct CubicSpline {
        int he = -1;        // "First" halfedge describing the canonical direction
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };
}   // namespace Curvenet