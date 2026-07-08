// dcurvenet_color.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <utility>

namespace DCurvenet {

    // Deformation
    struct vertColorData {
        // Runtime variables
        Eigen::Vector3d rgb;
    };

    struct heColorData {
        Eigen::Vector3d rgb;
    };

    struct faceDeformData {
    };
}   // namespace DCurvenet