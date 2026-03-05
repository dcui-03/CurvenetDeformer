// Utils.hpp
#pragma once

#include <Eigen/Core>
#include <glm/vec3.hpp>
#include <vector>

namespace Utils {
    // GLM::vec3 to Eigen::Vector3d converter
    Eigen::Vector3d glmToEigen(const glm::vec3 input);

    // Eigen::Vector3d to GLM::vec3 converter
    glm::vec3 eigenToGLM(const Eigen::Vector3d input);

    // Entire mesh conversion routine Eigen to GLM
    void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM);

    // Entire mesh conversion routine GLM to Eigen
    // Note that GLM is float, while Eigen prefers double
    // Do NOT convert back and forth, you will lose information
    void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM);

    // Copy positions and connectivity into a copied container
    void copyPositions(const std::vector<Eigen::Vector3d>& V_old, std::vector<Eigen::Vector3d>& V_new);

    void copyConnectivity(const std::vector<std::vector<int>>& T_old, std::vector<std::vector<int>>& T_new);


    // TODO: Fill in these functions
    // Project a vector onto a tangent plane, given the normal to the plane
    // Returns 1 if degenerate (shouldn't happen but we should handle it somehow?)
    int projectOntoTangentPlane(const Eigen::Vector3d& normal, Eigen::Vector3d& projection, bool normalize = true);

    // Given a set of projected vectors, order them CCW
    // TODO: How do we indicate degenerate vectors?
    Eigen::

} // namespace Utils