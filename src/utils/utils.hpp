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

// VECTOR/PROJECTION HELPERS

// Project a vector onto a tangent plane, given the normal to the plane
// Returns -1 if degenerate (shouldn't happen but we should handle it)
double projectVectorOntoTangentPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& vec, Eigen::Vector3d& proj, double scale = 1.0);

// Projects a point onto the tangent plane of a normal given a center 
Eigen::Vector3d projectPointOntoPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& center, const Eigen::Vector3d& p);

// Build basis (t1, t2) for a plane given a normal
void buildPlaneBasis(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2);

// Given a point on a plane basis and the plane basis, convert to 2D planar point
Eigen::Vector2d convertTo2D(const Eigen::Vector3d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2);

// Given a 2D planar point and the plane basis, revert to its 3D counterapart
Eigen::Vector3d revertTo3D(const Eigen::Vector2d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2);

// Check if a 2D point is in a 2D polygon
// To do this, we do raycasting to the segment
bool pointInPolygon2D(const Eigen::Vector2d& p, const std::vector<Eigen::Vector2d>& poly);

// Get the closest point on a segment in 2D and 3D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
// TODO: Can we combine the 2D and 3D cases using VectorXd?
Eigen::Vector2d closestPointOnSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, bool clip = true);

Eigen::Vector3d closestPointOnSegment3D(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, bool clip = true);
// Given a set of projected vectors, order them CCW
// TODO: How do we indicate degenerate vectors?

}  // namespace Utils
