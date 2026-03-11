#include "dvert.hpp"

#include <algorithm>

namespace DCurvenet {

dvert::dvert(const Eigen::Vector3d& position, int controlIndex, int parentSpline, double tValue)
    : pos(position), control(controlIndex), t(std::clamp(tValue, 0.0, 1.0)), parent_spline(parentSpline) {}

void dvert::addDSplineIdx(int dsplineIdx) { dspline_idxs.push_back(dsplineIdx); }

void dvert::setDSplineIdxs(const std::vector<int>& dsplineIdxs) { dspline_idxs = dsplineIdxs; }

void dvert::addSegmentIdx(int segmentIdx) { segment_idxs.push_back(segmentIdx); }

int dvert::numOutgoing() const {
    if (isControl()) {
        return static_cast<int>(dspline_idxs.size());
    }
    return static_cast<int>(segment_idxs.size());
}

}  // namespace DCurvenet
