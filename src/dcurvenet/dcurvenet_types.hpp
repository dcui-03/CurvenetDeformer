// dcurvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <utility>

// File with basic structs used by discrete curvenet class

namespace DCurvenet {

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();    // Init to zero, since most vertices will not receive an initial normal
        std::vector<int> adjHE;    // For control vertices, stores CCW outgoing HE's. For all others, stores a single outgoing halfedge
        bool active = true;     // For safety, say if the component is active (ignore for now)
        
        // cn_idx is redundant due to initialization, but better to be safe
        int cn_idx = -1;    // curvenet index if coincident with a control vertex
        int cn_type = -1;   // curvenet vertex type if coincident with a control vertex

        // Runtime variables
        Eigen::Vector3d new_pos;
    };

    // NOTE: A halfedge's next/prev can be -1 if this is the end of a spline
    // BUT if entering a high valence vertex, its next should be the CCW outgoing halfedge (and vice versa for prev)
    struct HalfEdge {
        int dest = -1;
        int twin = -1;
        int next = -1;
        int prev = -1;
        int edge = -1;
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Scaled Local Frame
        bool sign;  // true for positive, false for negative
        // Note, we can save these each separately for easy access.
        // Paper provides an easy method for computing def grad using components rather than matrices
        Eigen::Vector3d rest_tangent;
        Eigen::Vector3d rest_binormal;
        Eigen::Vector3d rest_normal;
        double rest_l, rest_w, rest_h; // length, width, and height

        // Runtime info: Altered frame needed for def grad computation
        Eigen::Vector3d tangent;
        Eigen::Vector3d binormal;
        Eigen::Vector3d normal;
        double l, w, h;
    };

    // Edge in a curve
    struct Edge {
        int he = -1;
        int curve = -1; // Index to curve
        bool active = true;
    };

    // Discrete curve
    struct Curve {
        int he_start = -1;        // Starting outgoing halfedge
        int he_end = -1;        // Ending outgoing halfedge
        int start = -1;    // One of the endpoint vertices
        int end = -1;      // The other endpoint vertex
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // NOTE: Not sure if corner attributes need to be kept as permanent fixtures. Maybe move to a temp variable instead?
        //       Right now, at runtime I'm just letting the new versions overwrite the old ones
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