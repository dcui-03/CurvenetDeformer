#pragma once

#include <Eigen/Core>

#include <vector>

namespace DCurvenet {

struct SegmentData {
    int dspline_idx = -1;
    int segment_idx = -1;
    Eigen::Vector3d t_out = Eigen::Vector3d::Zero();
    double l = 0.0;
    Eigen::Vector3d n_plus = Eigen::Vector3d::Zero();
    Eigen::Vector3d n_minus = Eigen::Vector3d::Zero();
    double w_plus = 0.0;
    double w_minus = 0.0;
};

class dvert {
public:
    dvert() = default;
    dvert(const Eigen::Vector3d& position, int controlIndex = -1, int parentSpline = -1, double tValue = 0.0);

    const Eigen::Vector3d& position() const { return pos; }
    int controlIndex() const { return control; }
    int parentSpline() const { return parent_spline; }
    double parameterT() const { return t; }
    bool isControl() const { return control >= 0; }

    void addDSplineIdx(int dsplineIdx);
    void setDSplineIdxs(const std::vector<int>& dsplineIdxs);
    const std::vector<int>& dsplineIdxs() const { return dspline_idxs; }

    void setCornerNormals(const std::vector<Eigen::Vector3d>& normals) { corner_normals = normals; }
    const std::vector<Eigen::Vector3d>& cornerNormals() const { return corner_normals; }
    void setSegmentData(const std::vector<SegmentData>& data) { segment_data = data; }
    const std::vector<SegmentData>& getSegmentData() const { return segment_data; }

    void addSegmentIdx(int segmentIdx);
    const std::vector<int>& segmentIdxs() const { return segment_idxs; }
    int numOutgoing() const;

private:
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    int control = -1;
    std::vector<int> dspline_idxs;
    std::vector<Eigen::Vector3d> corner_normals;
    std::vector<SegmentData> segment_data;
    std::vector<int> segment_idxs;
    double t = 0.0;
    int parent_spline = -1;
};

}  // namespace DCurvenet
