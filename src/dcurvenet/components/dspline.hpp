#pragma once

#include <vector>

namespace DCurvenet {

class dspline {
public:
    dspline() = default;
    dspline(int splineIdx, int startDvert, int endDvert);

    void setSegments(std::vector<int> segmentIdxs);
    void setInteriorSamples(int interiorSampleCount) { n_interior = interiorSampleCount; }

    const std::vector<int>& segmentIdxs() const { return segments; }
    int startDvert() const { return start_dvert; }
    int endDvert() const { return end_dvert; }
    int parentSpline() const { return parent_spline; }
    int interiorSamples() const { return n_interior; }
    int firstSegment() const { return first_segment; }
    int lastSegment() const { return last_segment; }

private:
    std::vector<int> segments;
    int first_segment = -1;
    int last_segment = -1;
    int start_dvert = -1;
    int end_dvert = -1;
    int n_interior = 0;
    int parent_spline = -1;
};

}  // namespace DCurvenet
