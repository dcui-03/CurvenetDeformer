// utils.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"

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
Eigen::Vector3d anyUnitTangent(const Eigen::Vector3d& normal);
Eigen::Vector3d anyPerpendicularUnit(const Eigen::Vector3d& tangent);
Eigen::Vector3d projectAndNormalizeToTangentPlane(const Eigen::Vector3d& normal,
                                                  const Eigen::Vector3d& tangent);
double computeMeanMeshEdgeLength(const std::vector<Eigen::Vector3d>& verts,
                                 const std::vector<std::vector<int>>& faces);

// Generic helpers used by main visualization path
void loadObjMesh(const std::string& path,
                 std::vector<Eigen::Vector3d>& vertices,
                 std::vector<std::vector<int>>& faces);

void buildPolyscopeCurveNetwork(const Curvenet::curvenet& cn,
                                std::vector<std::array<double, 3>>& points,
                                std::vector<std::array<std::size_t, 2>>& edges,
                                std::size_t samplesPerSpline);

void buildPolyscopeDiscreteCurveNetwork(const DCurvenet::dcurvenet& dcn,
                                        std::vector<std::array<double, 3>>& points,
                                        std::vector<std::array<std::size_t, 2>>& edges);

void buildPolyscopeControlCornerNormals(const DCurvenet::dcurvenet& dcn,
                                        std::vector<std::array<double, 3>>& origins,
                                        std::vector<std::array<double, 3>>& vectors);

std::vector<std::array<double, 3>> buildControlCloud(const Curvenet::curvenet& cn);

// TODO: Fill in when needed by the deformation stage.
// std::vector<int> sortVectorsCCW(...);

}  // namespace Utils
