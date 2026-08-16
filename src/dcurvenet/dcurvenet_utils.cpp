#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include <algorithm>

namespace DCurvenet {

    bool dcurvenet::isPositiveHalfedge(int he) const {
        return E[HE[he].edge].he == he;
    }

}   // namespace DCurvenet