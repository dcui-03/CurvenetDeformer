#include "pdcurvenet.hpp"

#include "igl/point_mesh_squared_distance.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace DCurvenet {

namespace {

std::uint64_t edgeKey(int a, int b) {
    const auto lo = static_cast<std::uint32_t>(std::min(a, b));
    const auto hi = static_cast<std::uint32_t>(std::max(a, b));
    return (static_cast<std::uint64_t>(lo) << 32U) | static_cast<std::uint64_t>(hi);
}

Eigen::Vector3d barycentricCoordinates(const Eigen::Vector3d& p,
                                      const Eigen::Vector3d& a,
                                      const Eigen::Vector3d& b,
                                      const Eigen::Vector3d& c) {
    const Eigen::Vector3d v0 = b - a;
    const Eigen::Vector3d v1 = c - a;
    const Eigen::Vector3d v2 = p - a;

    const double d00 = v0.dot(v0);
    const double d01 = v0.dot(v1);
    const double d11 = v1.dot(v1);
    const double d20 = v2.dot(v0);
    const double d21 = v2.dot(v1);
    const double denom = d00 * d11 - d01 * d01;

    if (std::abs(denom) <= std::numeric_limits<double>::epsilon()) {
        return Eigen::Vector3d(1.0, 0.0, 0.0);
    }

    const double v = (d11 * d20 - d01 * d21) / denom;
    const double w = (d00 * d21 - d01 * d20) / denom;
    const double u = 1.0 - v - w;
    return Eigen::Vector3d(u, v, w);
}

double pointSegmentSquaredDistance(const Eigen::Vector3d& p,
                                   const Eigen::Vector3d& a,
                                   const Eigen::Vector3d& b,
                                   double& t,
                                   Eigen::Vector3d& closest) {
    const Eigen::Vector3d ab = b - a;
    const double denom = ab.squaredNorm();
    if (denom <= std::numeric_limits<double>::epsilon()) {
        t = 0.0;
        closest = a;
        return (p - a).squaredNorm();
    }

    t = (p - a).dot(ab) / denom;
    t = std::clamp(t, 0.0, 1.0);
    closest = a + t * ab;
    return (p - closest).squaredNorm();
}

}  // namespace

pdcurvenet::pdcurvenet(const dcurvenet& parentDC,
                       const std::vector<Eigen::Vector3d>& meshV,
                       const std::vector<std::vector<int>>& meshF,
                       const Eigen::MatrixXi& triF,
                       const Eigen::VectorXi& triToFace,
                       double eps)
    : dcurvenet(parentDC),
      mesh_vertices_(meshV),
      mesh_faces_(meshF),
      tri_f_matrix_(triF) {
    tri_to_mesh_face_.assign(triToFace.data(), triToFace.data() + triToFace.size());

    triangle_faces_.resize(static_cast<std::size_t>(tri_f_matrix_.rows()));
    for (int i = 0; i < tri_f_matrix_.rows(); ++i) {
        triangle_faces_[static_cast<std::size_t>(i)] = {
            tri_f_matrix_(i, 0), tri_f_matrix_(i, 1), tri_f_matrix_(i, 2)};
    }

    mesh_v_matrix_.resize(static_cast<int>(mesh_vertices_.size()), 3);
    for (std::size_t i = 0; i < mesh_vertices_.size(); ++i) {
        mesh_v_matrix_.row(static_cast<int>(i)) = mesh_vertices_[i];
    }

    buildMeshEdgeTable();
    ComputePDC(eps);
}

void pdcurvenet::ComputePDC(double eps) {
    if (mesh_vertices_.empty()) {
        throw std::invalid_argument("pdcurvenet requires a non-empty mesh vertex list");
    }
    if (mesh_faces_.empty() || triangle_faces_.empty()) {
        throw std::invalid_argument("pdcurvenet requires a mesh with at least one valid face");
    }

    snap_epsilon_ = (eps > 0.0) ? eps : computeDefaultSnapEpsilon();
    projectToSurface();
    snapToNearbyElements();
}


void pdcurvenet::buildMeshEdgeTable() {
    mesh_edges_.clear();
    std::unordered_map<std::uint64_t, int> edge_to_index;

    for (std::size_t face_idx = 0; face_idx < mesh_faces_.size(); ++face_idx) {
        const auto& face = mesh_faces_[face_idx];
        if (face.size() < 2) {
            continue;
        }

        for (std::size_t i = 0; i < face.size(); ++i) {
            const int a = face[i];
            const int b = face[(i + 1) % face.size()];
            const std::uint64_t key = edgeKey(a, b);

            auto [it, inserted] = edge_to_index.emplace(key, static_cast<int>(mesh_edges_.size()));
            if (inserted) {
                mesh_edges_.push_back({a, b, {}});
            }
            mesh_edges_[static_cast<std::size_t>(it->second)].incident_faces.push_back(static_cast<int>(face_idx));
        }
    }
}

void pdcurvenet::projectToSurface() {
    const auto& dverts = verts();
    projected_samples_.assign(dverts.size(), {});

    Eigen::MatrixXd queries(static_cast<int>(dverts.size()), 3);
    for (std::size_t i = 0; i < dverts.size(); ++i) {
        queries.row(static_cast<int>(i)) = dverts[i].position();
    }

    Eigen::VectorXd sqr_d;
    Eigen::VectorXi tri_ids;
    Eigen::MatrixXd closest_points;
    igl::point_mesh_squared_distance(queries, mesh_v_matrix_, tri_f_matrix_, sqr_d, tri_ids, closest_points);

    for (std::size_t i = 0; i < dverts.size(); ++i) {
        ProjectedSample sample;
        sample.dvert_index = static_cast<int>(i);
        sample.source_position = dverts[i].position();
        sample.projected_position = closest_points.row(static_cast<int>(i));
        sample.squared_distance = sqr_d(static_cast<int>(i));
        sample.triangle_index = tri_ids(static_cast<int>(i));

        if (sample.triangle_index < 0 ||
            static_cast<std::size_t>(sample.triangle_index) >= triangle_faces_.size()) {
            throw std::runtime_error("libigl returned an invalid closest-triangle index");
        }

        sample.mesh_face = tri_to_mesh_face_[static_cast<std::size_t>(sample.triangle_index)];
        sample.triangle_vertices = triangle_faces_[static_cast<std::size_t>(sample.triangle_index)];

        const Eigen::Vector3d& a = mesh_vertices_[static_cast<std::size_t>(sample.triangle_vertices[0])];
        const Eigen::Vector3d& b = mesh_vertices_[static_cast<std::size_t>(sample.triangle_vertices[1])];
        const Eigen::Vector3d& c = mesh_vertices_[static_cast<std::size_t>(sample.triangle_vertices[2])];
        sample.triangle_barycentric = barycentricCoordinates(sample.projected_position, a, b, c);

        projected_samples_[i] = sample;
    }
}

void pdcurvenet::snapToNearbyElements() {
    const double eps2 = snap_epsilon_ * snap_epsilon_;

    std::unordered_map<std::uint64_t, int> edge_to_index;
    edge_to_index.reserve(mesh_edges_.size());
    for (std::size_t edge_idx = 0; edge_idx < mesh_edges_.size(); ++edge_idx) {
        const auto& edge = mesh_edges_[edge_idx];
        edge_to_index.emplace(edgeKey(edge.v0, edge.v1), static_cast<int>(edge_idx));
    }

    for (auto& sample : projected_samples_) {
        if (sample.mesh_face < 0 || static_cast<std::size_t>(sample.mesh_face) >= mesh_faces_.size()) {
            continue;
        }

        const auto& face = mesh_faces_[static_cast<std::size_t>(sample.mesh_face)];

        double best_vertex_d2 = std::numeric_limits<double>::infinity();
        int best_vertex = -1;
        for (int vid : face) {
            const double d2 = (sample.projected_position - mesh_vertices_[static_cast<std::size_t>(vid)]).squaredNorm();
            if (d2 < best_vertex_d2) {
                best_vertex_d2 = d2;
                best_vertex = vid;
            }
        }

        if (best_vertex >= 0 && best_vertex_d2 <= eps2) {
            sample.attachment = MeshAttachmentType::Vertex;
            sample.mesh_vertex = best_vertex;
            sample.mesh_edge = -1;
            sample.projected_position = mesh_vertices_[static_cast<std::size_t>(best_vertex)];
            sample.edge_t = 0.0;
            continue;
        }

        double best_edge_d2 = std::numeric_limits<double>::infinity();
        int best_edge = -1;
        double best_t = 0.0;
        Eigen::Vector3d best_point = sample.projected_position;

        for (std::size_t i = 0; i < face.size(); ++i) {
            const int a = face[i];
            const int b = face[(i + 1) % face.size()];

            double t = 0.0;
            Eigen::Vector3d closest = Eigen::Vector3d::Zero();
            const double d2 = pointSegmentSquaredDistance(
                sample.projected_position,
                mesh_vertices_[static_cast<std::size_t>(a)],
                mesh_vertices_[static_cast<std::size_t>(b)],
                t,
                closest);
            if (d2 < best_edge_d2) {
                best_edge_d2 = d2;
                best_t = t;
                best_point = closest;

                const auto it = edge_to_index.find(edgeKey(a, b));
                best_edge = (it == edge_to_index.end()) ? -1 : it->second;
            }
        }

        if (best_edge >= 0 && best_edge_d2 <= eps2) {
            sample.attachment = MeshAttachmentType::Edge;
            sample.mesh_vertex = -1;
            sample.mesh_edge = best_edge;
            sample.projected_position = best_point;
            sample.edge_t = best_t;
            continue;
        }

        sample.attachment = MeshAttachmentType::Face;
        sample.mesh_vertex = -1;
        sample.mesh_edge = -1;
        sample.edge_t = 0.0;
    }
}

double pdcurvenet::computeDefaultSnapEpsilon() const {
    if (mesh_vertices_.empty()) {
        return 1e-6;
    }

    Eigen::Vector3d bbox_min = mesh_vertices_.front();
    Eigen::Vector3d bbox_max = mesh_vertices_.front();
    for (const auto& v : mesh_vertices_) {
        bbox_min = bbox_min.cwiseMin(v);
        bbox_max = bbox_max.cwiseMax(v);
    }

    const double diagonal = (bbox_max - bbox_min).norm();
    return std::max(1e-6, 1e-5 * diagonal);
}

}  // namespace DCurvenet
