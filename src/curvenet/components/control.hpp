// control.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

namespace Curvenet {


// Control Points
class control {
    public:
        control(Eigen::Vector3d controlP, Eigen::Vector3d normal);
        
        // TODO: Needs getters so that the UI can ask for internals
        Eigen::Vector3d getPos();
        void setPos(Eigen::Vector3d new_pos);

        Eigen::Vector3d getNormal();

        std::vector<int> getSplineIdxs();
        void setSplineIdxs(std::vector<int> reordering);
        void addSplineIdx(int s_idx);
    protected:
        // No class inheritance
    private:
        // Number of outgoing splines from this control
        int num_outgoing();

        // Store attributes
        Eigen::Vector3d pos;
        Eigen::Vector3d normal; // Associated surface normal
        std::vector<int> spline_idxs;   // pointers to an ordering of outgoing splines CCW
        // Do we need to explicitly store tangents? --> No, probably not.
};

// Tangent Points
class tangent {
    public:
        tangent(Eigen::Vector3d tangentP, int control_idx, int spline_idx);
        Eigen::Vector3d getPos();
        void setPos(Eigen::Vector3d new_pos);
    protected:
        // No class inheritance
    private:
        // Store attributes
        Eigen::Vector3d pos;
        int parent_control;   // pointer to associated control
        int spline_idx;       // pointer to associated spline
};

}   // namespace Curvenet