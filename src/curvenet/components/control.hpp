// control.hpp
#pragma once

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace Curvenet {

// Control Points
class control {
public:

    control(Eigen::Vector3d controlP, Eigen::Vector3d normal);

    // TODO: Needs getters so that the UI can ask for internals
    const Eigen::Vector3d& getPosition() const { return pos; }
    const Eigen::Vector3d& getNormal() const { return normal; }

    void setPosition(const Eigen::Vector3d& p) { pos = p; }
    void addSplineIdx(int splineIdx) { spline_idxs.push_back(splineIdx); }
    void setSplineIdxs(const std::vector<int>& idxs) { spline_idxs = idxs; }

    // Number of outgoing splines from this control
    int num_outgoing() const;

    const std::vector<int>& getSplineIdxs() const { return spline_idxs; }

protected:
    // No class inheritance
private:
    // Store attributes
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    std::vector<int> spline_idxs;  // pointers to an ordering of outgoing splines CCW
};

// Tangent Points
class tangent {
public:

    tangent(Eigen::Vector3d controlP, int parentControl = -1);

    const Eigen::Vector3d& getPosition() const { return pos; }
    int getParentControl() const { return parent_control; }
    int getSplineIdx() const { return spline_idx; }

    void setPosition(const Eigen::Vector3d& p) { pos = p; }
    void setParentControl(int idx) { parent_control = idx; }
    void setSplineIdx(int idx) { spline_idx = idx; }

protected:
    // No class inheritance
private:
    // Store attributes
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    int parent_control = -1;  // pointer to associated control
    int spline_idx = -1;      // pointer to associated spline
};

}  // namespace Curvenet
