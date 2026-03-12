#include "control.hpp"

#include <Eigen/Core>
#include <vector>

namespace Curvenet {
    // CONTROLS
    control::control(Eigen::Vector3d controlP, Eigen::Vector3d normal):
        pos(controlP), normal(normal) {
    }

    Eigen::Vector3d control::getPos() {
        return pos;
    }

    void control::setPos(Eigen::Vector3d new_pos) {
        pos = new_pos;
    }

    Eigen::Vector3d control::getNormal() {
        return normal;
    }

    int control::num_outgoing() {
        return spline_idxs.size();
    }

    // Get a copy of the outgoing spline list
    std::vector<int> control::getSplineIdxs() {
        return spline_idxs;
    }

    void control::setSplineIdxs(std::vector<int> reordering) {
        spline_idxs = reordering;
        return;
    }

    // Add a spline to the control's list
    void control::addSplineIdx(int s_idx) {
        spline_idxs.push_back(s_idx);
        return;
    }

    // TANGENTS
    tangent::tangent(Eigen::Vector3d tangentP, int control_idx, int spline_idx):
        pos(tangentP), parent_control(control_idx), spline_idx(spline_idx) {
    }

    Eigen::Vector3d tangent::getPos() {
        return pos;
    }

    void tangent::setPos(Eigen::Vector3d new_pos) {
        pos = new_pos;
    }

}   // namespace Curvenet