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

        // Other indices for cut mesh
        int dCN_idx = -1;
        bool SG_vert = false;
    };

    struct HalfEdge {
        int dest = -1;
        int twin = -1;
        int next = -1;
        int prev = -1;
        int edge = -1;           // Index to the edge
        int face = -1;           // Index to the face
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Other indices for cut mesh
        int dCN_idx = -1;   // -1 if not connected
        bool dCN_sign = true;  // true for positive, false for negative
    };

    // Only store one halfedge for each edge
    // No real need to store the normal
    struct Edge {
        int he = -1;
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    struct Face {
        int he = -1;
        std::vector<int> verts; // Face vertices. NOTE: These are NOT necessarily in order!
        Eigen::Vector3d n = Eigen::Vector3d::Zero();  // Face normal
        double fArea = 0.0;   // Face area. To be computed only when necessary
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };
}   // namespace Mesh