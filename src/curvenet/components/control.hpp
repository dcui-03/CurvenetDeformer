// control.hpp
#pragma once

#include "spline.hpp"
#include <Eigen/Core>
#include <vector>

namespace Curvenet {


// Control Points
class control {
    public:
        control(Eigen::Vector3d controlP, Eigen::Vector3d);

        // TODO: Needs getters so that the UI can ask for internals
    protected:
        // No class inheritance
    private:
        // Number of outgoing splines from this control
        int num_outgoing();

        // Store attributes
        Eigen::Vector3d pos;
        Eigen::Vector3d normal;
        std::vector<int> spline_idxs;   // pointers to an ordering of outgoing splines CCW
        
};

// Tangent Points
class tangent {
    public:
        tangent(Eigen::Vector3d controlP);
    protected:
        // No class inheritance
    private:
        // Store attributes
        Eigen::Vector3d pos;
        int parent_control;   // pointer to associated control
        int spline_idx;       // pointer to associated spline
}

}   // namespace Curvenet