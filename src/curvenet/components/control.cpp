#include "control.hpp"

#include <utility>

namespace Curvenet {

control::control(Eigen::Vector3d controlP, Eigen::Vector3d controlNormal)
    : pos(std::move(controlP)), normal(std::move(controlNormal)) {}

int control::num_outgoing() const { return static_cast<int>(spline_idxs.size()); }

tangent::tangent(Eigen::Vector3d controlP, int parentControl)
    : pos(std::move(controlP)), parent_control(parentControl) {}

}  // namespace Curvenet
