// control.hpp
#pragma once

#include <Eigen/Core>

#include <vector>

namespace Curvenet {

class control {
public:
    control(Eigen::Vector3d controlP, Eigen::Vector3d normal);

    const Eigen::Vector3d& getPosition() const { return pos; }
    const Eigen::Vector3d& getPos() const { return pos; }
    const Eigen::Vector3d& getNormal() const { return normal; }

    void setPosition(const Eigen::Vector3d& p) { pos = p; }
    void setPos(const Eigen::Vector3d& p) { pos = p; }

    void addSplineIdx(int splineIdx) { spline_idxs.push_back(splineIdx); }
    void setSplineIdxs(const std::vector<int>& idxs) { spline_idxs = idxs; }
    const std::vector<int>& getSplineIdxs() const { return spline_idxs; }

    int num_outgoing() const;

private:
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    std::vector<int> spline_idxs;
};

class tangent {
public:
    tangent(Eigen::Vector3d tangentP, int parentControl = -1, int splineIdx = -1);

    const Eigen::Vector3d& getPosition() const { return pos; }
    const Eigen::Vector3d& getPos() const { return pos; }
    int getParentControl() const { return parent_control; }
    int getSplineIdx() const { return spline_idx; }

    void setPosition(const Eigen::Vector3d& p) { pos = p; }
    void setPos(const Eigen::Vector3d& p) { pos = p; }
    void setParentControl(int idx) { parent_control = idx; }
    void setSplineIdx(int idx) { spline_idx = idx; }

private:
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    int parent_control = -1;
    int spline_idx = -1;
};

}  // namespace Curvenet
