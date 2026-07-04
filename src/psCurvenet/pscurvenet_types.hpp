// pscurvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace psCurvenet {
    // NOTE: These use pointers, which should make deletion much easier
    struct Control {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        std::vector<HalfEdge*> adjHE; // Adjacent halfedges in no particular order
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // NOTE: Halfedge iteration indices (next, prev) prioritize easy iteration over their associated controls
    //       As a result, they may not be ideal for standard halfedge traversal tasks.
    struct HalfEdge {
        Control* origin;  // To control
        Spline* s;     // associated spline
        // Iterators: Note that if we are at an endpoint, leave entry as -1
        HalfEdge* twin;
        Eigen::Vector3d tan;   // Tangent vector (defined in global coordinates, NOT relative to control)

        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Splines for easy iteration
    struct Spline {
        Control* start;
        Control* end;
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

}   // namespace psCurvenet