// utils.hpp
#pragma once

#include "curvenet/curvenet.hpp"

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace Utils {

#if 0
// Legacy conversion helpers preserved from original stub.
// Disabled because GLM is not a required dependency in the current build.
Eigen::Vector3d glmToEigen(const glm::vec3 input);
glm::vec3 eigenToGLM(const Eigen::Vector3d input);
void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM);
void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM);
#endif

// Copy positions and connectivity into a copied container
void copyPositions(const std::vector<Eigen::Vector3d>& oldV, std::vector<Eigen::Vector3d>& newV);
void copyConnectivity(const std::vector<std::vector<int>>& oldT, std::vector<std::vector<int>>& newT);

// Generic helpers used by main visualization path
void loadObjMesh(const std::string& path,
                 std::vector<Eigen::Vector3d>& vertices,
                 std::vector<std::vector<int>>& faces);

void buildPolyscopeCurveNetwork(const Curvenet::curvenet& cn,
                                std::vector<std::array<double, 3>>& points,
                                std::vector<std::array<std::size_t, 2>>& edges,
                                std::size_t samplesPerSpline);

std::vector<std::array<double, 3>> buildControlCloud(const Curvenet::curvenet& cn);

// TODO: Fill in when needed by the deformation stage.
// int projectOntoTangentPlane(const Eigen::Vector3d& normal, Eigen::Vector3d& projection, bool normalize = true);
// std::vector<int> sortVectorsCCW(...);

}  // namespace Utils
