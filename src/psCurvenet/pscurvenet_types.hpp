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
        Eigen::Vector3d bn = Eigen::Vector3d::Zero();    // A "Canonical" basis for the control to enable in-plane tangent rotations
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Splines for easy iteration
    struct Spline {
        int start;
        int end;
        Eigen::Vector3d t0;
        Eigen::Vector3d t1;
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

}   // namespace psCurvenet