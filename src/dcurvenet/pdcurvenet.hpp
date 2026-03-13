// pdcurvenet.hpp
#pragma once

#include "dcurvenet.hpp"

#include <Eigen/Core>

#include <array>
#include <limits>
#include <vector>

namespace DCurvenet {

enum class MeshAttachmentType {
    Vertex,
    Edge,
    Face
};

struct MeshEdgeInfo {
    int v0 = -1;
    int v1 = -1;
    std::vector<int> incident_faces;
};

struct ProjectedSample {
    int dvert_index = -1;
    Eigen::Vector3d source_position = Eigen::Vector3d::Zero();
    Eigen::Vector3d projected_position = Eigen::Vector3d::Zero();
    double squared_distance = std::numeric_limits<double>::infinity();

    MeshAttachmentType attachment = MeshAttachmentType::Face;
    int mesh_vertex = -1;
    int mesh_edge = -1;
    int mesh_face = -1;

    int triangle_index = -1;
    std::array<int, 3> triangle_vertices{{-1, -1, -1}};
    Eigen::Vector3d triangle_barycentric = Eigen::Vector3d::Zero();

    double edge_t = 0.0;
};

class pdcurvenet : public dcurvenet {
public:
    pdcurvenet() = default;
    pdcurvenet(const dcurvenet& parentDC,
               const std::vector<Eigen::Vector3d>& meshV,
               const std::vector<std::vector<int>>& meshF,
               double eps = -1.0);

    void ComputePDC(double eps = -1.0);

    const std::vector<Eigen::Vector3d>& meshVertices() const { return mesh_vertices_; }
    const std::vector<std::vector<int>>& meshFaces() const { return mesh_faces_; }
    const std::vector<std::array<int, 3>>& triangleFaces() const { return triangle_faces_; }
    const std::vector<int>& triangleToMeshFace() const { return tri_to_mesh_face_; }
    const std::vector<MeshEdgeInfo>& meshEdges() const { return mesh_edges_; }
    const std::vector<ProjectedSample>& projectedSamples() const { return projected_samples_; }
    double snapEpsilon() const { return snap_epsilon_; }

private:
    void buildTriangulatedSurface();
    void buildMeshEdgeTable();
    void projectToSurface();
    void snapToNearbyElements();
    double computeDefaultSnapEpsilon() const;

    std::vector<Eigen::Vector3d> mesh_vertices_;
    std::vector<std::vector<int>> mesh_faces_;
    std::vector<std::array<int, 3>> triangle_faces_;
    std::vector<int> tri_to_mesh_face_;
    std::vector<MeshEdgeInfo> mesh_edges_;
    std::vector<ProjectedSample> projected_samples_;

    Eigen::MatrixXd mesh_v_matrix_;
    Eigen::MatrixXi tri_f_matrix_;

    double snap_epsilon_ = 0.0;
};

}  // namespace DCurvenet
