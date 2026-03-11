#pragma once

#include <Eigen/Core>

namespace DCurvenet {

struct SegmentSideData {
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    Eigen::Vector3d b = Eigen::Vector3d::Zero();
    double w = 0.0;
    double h = 0.0;
    Eigen::Matrix3d B = Eigen::Matrix3d::Identity();
    Eigen::Vector3d S = Eigen::Vector3d::Zero();
    Eigen::Matrix3d BS = Eigen::Matrix3d::Zero();
    bool valid = false;
};

class dsegment {
public:
    dsegment() = default;
    dsegment(int startIdx,
             int endIdx,
             int dsplineIdx,
             const Eigen::Vector3d& startPos,
             const Eigen::Vector3d& endPos);

    int startDvert() const { return start_dvert; }
    int endDvert() const { return end_dvert; }
    int parentDSpline() const { return parent_dspline; }

    const Eigen::Vector3d& direction() const { return direc; }
    double segmentLength() const { return length; }
    const SegmentSideData& plusSide() const { return plus_side; }
    const SegmentSideData& minusSide() const { return minus_side; }
    void setPlusSide(const SegmentSideData& side) { plus_side = side; }
    void setMinusSide(const SegmentSideData& side) { minus_side = side; }

private:
    Eigen::Vector3d direc = Eigen::Vector3d::Zero();
    double length = 0.0;

    int start_dvert = -1;
    int end_dvert = -1;
    int parent_dspline = -1;

    SegmentSideData plus_side;
    SegmentSideData minus_side;
};

}  // namespace DCurvenet
