// straightestGeodesic.hpp
#pragma once

#include "dcurvenet/pdcurvenet.hpp"
#include "dcurvenet/components/dvert.hpp"
#include "utils.hpp"
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

// Implementation of Polthier & Schmies [1997], Straightest Geodesics on Polyhedral Surfaces
namespace StraightestGeodesics {
    // Return an ordered list of dverts from the start point to the end point
    // Note: We can get a starting direction by projecting onto the local tangent plane (get projected element using dvert)
    // At each split, perform snapping and reset movement direction using vector; if snapping to vert, check adjacent edge directions first
    // TODO: include mesh as input
    std::vector<DCurvenet::dvert> computeStraightestGeodesic(const DCurvenet::dvert& start_dvert, const DCurvenet::dvert& end_dvert, bool snap_verts = true);

    // Find intersection of a ray with edges on a face given a starting vertex on a face, edge or vertex
    // Returns the t value along the CCW oriented edge that is intersected.
    // Note: if we want to get the t-value, we can't do ray-plane intersection
    // We need to rotate into x-y plane and then unrotate after computing intersection
    double rayEdgeIntersectionF(std::vector<Eigen::Vector3d>& f, int start_edge, Eigen::Vector3d& direction);

    double rayEdgeIntersectionE(std::vector<Eigen::Vector3d>& f, int start_edge, Eigen::Vector3d& direction, double start_t);

    double rayEdgeIntersectionV(std::vector<Eigen::Vector3d>& f, int start_vert, Eigen::Vector3d& direction);

    // Project face vertices onto the plane spanned by the face normal (given a fixed edge)
    std::vector<Eigen::Vector3d> projectFaceUsingNormal(std::vector<Eigen::Vector3d>& f, int edge_idx, const Eigen::Vector3d& fNormal);

    // Rotation matrix to get a projected face into a tangent plane by edge
    // Takes the normal we want to match with the face normal we currently have
    Eigen::Matrix3d rotateAroundEdge(const Eigen::Vector3d& rAxis, const Eigen::Vector3d& fNormal, const Eigen::Vector3d targetNormal);



}   // namespace StraightestGeodesics