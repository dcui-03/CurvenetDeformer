// jsonUtils.hpp
#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

namespace JSONUtils {

struct CurvenetInput {
    std::vector<Eigen::Vector3d> controlP;
    std::vector<Eigen::Vector3d> surfaceN;
    std::vector<std::vector<int>> curveC;  // [start, tangent1, tangent2, end]
};

CurvenetInput loadBezierCurvenetInput(const std::string& jsonPath, double dedupTolerance = 1e-6);

}  // namespace JSONUtils

