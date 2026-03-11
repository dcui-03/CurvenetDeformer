#include "dspline.hpp"

#include <utility>

namespace DCurvenet {

dspline::dspline(int splineIdx, int startDvert, int endDvert)
    : start_dvert(startDvert), end_dvert(endDvert), parent_spline(splineIdx) {}

void dspline::setSegments(std::vector<int> segmentIdxs) {
    segments = std::move(segmentIdxs);
    if (segments.empty()) {
        first_segment = -1;
        last_segment = -1;
        return;
    }

    first_segment = segments.front();
    last_segment = segments.back();
}

}  // namespace DCurvenet
