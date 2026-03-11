#include "dsegment.hpp"

#include <limits>

namespace DCurvenet {

dsegment::dsegment(int startIdx,
                   int endIdx,
                   int dsplineIdx,
                   const Eigen::Vector3d& startPos,
                   const Eigen::Vector3d& endPos)
    : start_dvert(startIdx), end_dvert(endIdx), parent_dspline(dsplineIdx) {
    const Eigen::Vector3d diff = endPos - startPos;
    length = diff.norm();
    if (length > std::numeric_limits<double>::epsilon()) {
        direc = diff / length;
    }
}

}  // namespace DCurvenet
