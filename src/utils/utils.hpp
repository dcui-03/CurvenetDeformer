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

    // SORTING
    void sortAscending_IDsUsingValues(std::vector<int> idxs, std::vector<double> values);

    // Insert at index between a pair of indices in a list
    bool insertIdxBetweenPair(std::vector<int>& idxList, int a, int b, int new_idx);
    

    // VECTOR/PROJECTION HELPERS
    
    // Computes the closest point to a triangle
    Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d> triVerts, const Eigen::Vector3d p);

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

    bool raycastToSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& direc, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, double& t, double& u, bool clip = true);

    // Get the closest point on a segment in 2D and 3D, where the endpoints are defined
    // To do this, project onto parameterized segment and snap t to [0, 1]
    // TODO: Can we combine the 2D and 3D cases using VectorXd?
    Eigen::Vector2d closestPointOnSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, bool clip = true);

    Eigen::Vector3d closestPointOnSegment3D(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, bool clip = true);

    // Given points on a face, plus a ray, compute the intersection of the ray with the face, if one exists
    bool computeFaceIntersection(const std::vector<Eigen::Vector3d>& fVerts, const Eigen::Vector3d& fNormal,
                             const Eigen::Vector3d& start, const Eigen::Vector3d& direc,
                             double& t, int& el_type, int& local_idx, double& u, double& theta, double tol = 1e-6);
    bool computeFaceIntersectionTarget(const std::vector<Eigen::Vector3d>& fVerts, const Eigen::Vector3d& fNormal,
                             const Eigen::Vector3d& start, const Eigen::Vector3d& target,
                             double& t, int& el_type, int& local_idx, double& u, double& theta, double tol = 1e-6);

    // MEAN VALUE COORDINATES
    // Helpers
    double vectorAngle(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3);
    double computeSign(const double& value);

    // Returns the weights only
    // Follows method of Fuda and Hormann [2024]
    void meanValueCoordinates(const Eigen::Vector2d& target, const std::vector<Eigen::Vector2d>& cage, Eigen::VectorXd& weights);

} // namespace Utils