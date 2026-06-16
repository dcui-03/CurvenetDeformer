// dcurvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

// File with basic structs used by discrete curvenet class

namespace DCurvenet {

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        std::vector<int> adjHE;    // CCW ordered list of outgoing HE's
        bool active = true;     // For safety, say if the component is active (ignore for now)
        
        // Other info
        int cn_idx = -1;    // curvenet index if coincident with a control vertex
    };

    // NOTE: A halfedge's next/prev can be -1 if this is the end of a spline
    // BUT if entering a high valence vertex, its next should be the CCW outgoing halfedge (and vice versa for prev)
    struct HalfEdge {
        int dest = -1;
        int twin = -1;
        int next = -1;
        int prev = -1;
        int edge = -1;           // Index to the edge
        int spline = -1;           // Index to the spline
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
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Curve (discrete spline) for easy iteration
    struct Curve {
        int he = -1;        // Starting outgoing halfedge
        int start_v = -1;    // One of the endpoint vertices
        int end_v = -1;      // The other endpoint vertex
        double length = 0.0;   // Discrete spline length
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Other info
        int cn_idx = -1;    // Curvenet index if parent is a spline
    };
}   // namespace DCurvenet