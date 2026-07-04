// mesh_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

// File with basic structs used by mesh class

namespace Mesh {

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        double vArea = 0.0;   // barycentric dual area. To be computed only when necessary
        int he = -1;   // one outgoing halfedge, or -1 if isolated
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Attributes for cut mesh
        int label = 0;  // {0 if original mesh vertex, 1 if projected CN vertex, 2 otherwise}
        int dCN_idx = -1;   // Corresponding dCN index for cut-mesh (if cut-vertex is associated with a dCN vert or HE)
        // Note that if label == 1, then dCN_idx is a vert index, and if label == 2, then dCN_idx is a HE index
        // Where we landed in the original mesh, if label == 1
        int mesh_elType = -1;
        int mesh_elIdx = -1;
        // Vector from the projected point on the rest mesh to the rest curvenet vert (if label == 1)
        Eigen::Vector3d projVector = Eigen::Vector3d::Zero();
    };

    struct HalfEdge {
        int dest = -1;

        int twin = -1;
        int next = -1;
        int prev = -1;

        int edge = -1;           // Index to the edge
        int face = -1;           // Index to the face

        bool active = true;     // For safety, say if the component is active (ignore for now)
        bool boundary = false;  // This is necessary for cutmesh face reinitialization

        // Other indices for cut mesh
        int dCN_idx = -1;   // -1 if not connected, dCN HE index otherwise
        // Computation needs
        // Deformation gradient eventually computed using Laplacian
        Eigen::Matrix3d defGrad = Eigen::Matrix3d::Identity();
    };

    // Only store one halfedge for each edge
    // No real need to store the normal
    struct Edge {
        int he = -1;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    struct Face {
        int he = -1;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();  // Face normal
        double fArea = 0.0;   // Face area. To be computed only when necessary
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };
}   // namespace Mesh