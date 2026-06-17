// dcurvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

// File with basic structs used by discrete curvenet class

namespace DCurvenet {

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        std::vector<int> adjHE;    // For control vertices, stores CCW outgoing HE's. For all others, stores a single outgoing halfedge
        bool active = true;     // For safety, say if the component is active (ignore for now)
        
        // cn_idx is redundant due to initialization, but better to be safe
        int cn_idx = -1;    // curvenet index if coincident with a control vertex
        int cn_type = -1;   // curvenet vertex type if coincident with a control vertex
        // Vector from the projected point on the rest mesh to the curvenet vert (computed during cut-mesh computation)
        Eigen::Vector3d projVector_pos = Eigen::Vector3d::Zero();
        Eigen::Vector3d projVector_neg = Eigen::Vector3d::Zero();
    };

    // NOTE: A halfedge's next/prev can be -1 if this is the end of a spline
    // BUT if entering a high valence vertex, its next should be the CCW outgoing halfedge (and vice versa for prev)
    struct HalfEdge {
        int dest = -1;
        int twin = -1;
        int next = -1;
        int prev = -1;
        int edge = -1;           // Index to the edge
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Scaled Local Frame
        bool sign;  // true for positive, false for negative
        Eigen::Vector3d tangent;
        Eigen::Vector3d binormal;
        Eigen::Vector3d normal;
        Eigen::Vector3d scale = {1.0, 1.0, 1.0}; // Scale each axis
    };

    // Only store one halfedge for each edge
    struct Edge {
        int he = -1;
        int curve = -1;
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Discrete curve
    struct Curve {
        int he = -1;        // Starting outgoing halfedge
        int start = -1;    // One of the endpoint vertices
        int end = -1;      // The other endpoint vertex
        double len = 0.0;   // length of the curve
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Corner normals
        std::pair<Eigen::Vector3d, Eigen::Vector3d> N_pos;  // first is start, second is end
        std::pair<Eigen::Vector3d, Eigen::Vector3d> N_neg;  // first is start, second is end
        // Corner widths
        std::pair<double, double> W_pos;    // first is start, second is end
        std::pair<double, double> W_neg;    // first is start, second is end

        // Other info
        int cn_idx = -1;    // Curvenet index if parent is a spline
    };
}   // namespace DCurvenet