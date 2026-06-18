#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <cmath>

namespace DCurvenet {

    Eigen::Matrix3d dcurvenet::computeHEDefGrad(int he,
                                                const Eigen::Vector3d& newT, 
                                                const Eigen::Vector3d& newB, 
                                                const Eigen::Vector3d& newN, 
                                                const double& newL, const double& newW, const double& newH) {
        const Eigen::Vector3d& T = HE[he].tangent;
        const Eigen::Vector3d& B = HE[he].binormal;
        const Eigen::Vector3d& N = HE[he].normal;
        const double& l = HE[he].l;
        const double& w = HE[he].w;
        const double& h = HE[he].h;
        Eigen::Matrix3d F = (newL / l) * (newT * T.transpose()) + 
                            (newW / w) * (newB * B.transpose()) + 
                            (newH / h) * (newN * N.transpose());
        return F;
    }

    // Accumulate rotation matrices, starting from an initial halfedge and tracing forward until we hit the goal vertex
    double dcurvenet::accumulateRotations(int start_he, int end_v, std::vector<Eigen::Matrix3d>& rots, std::vector<double> lens) {
        rots.clear();
        lens.clear();
        int he_curr = start_he;
        int he_prev = start_he;
        Eigen::Vector3d tan_curr = HE[he_curr].tangent;
        Eigen::Vector3d tan_prev = tan_curr;
        std::vector<double> lengths;
        std::vector<Eigen::Matrix3d> rots;
        Eigen::Matrix3d curr_rot = Eigen::Matrix3d::Identity();
        double curr_len = 0.0;
        // Traverse halfedges until we hit the end vertex
        do {
            curr_rot = Utils::computeRotation(tan_prev, tan_curr) * curr_rot;   // Left multiply to accumulate rotations
            rots.push_back(curr_rot);
            curr_len += HE[he_curr].l;
            lengths.push_back(curr_len);

            he_prev = he_curr;
            he_curr = HE[he_curr].next;
            tan_prev = tan_curr;
            tan_curr = HE[he_curr].tangent;
        } while ((HE[he_prev].dest != end_v) || (he_curr == -1));
        return curr_len;
    }

    // Total torsion
    double dcurvenet::computeTorsion(Eigen::Vector3d n_1, Eigen::Vector3d n_k, Eigen::Matrix3d Om_k, Eigen::Vector3d t_k) {
        Eigen::Vector3d twist = Om_k * n_1;
        return std::atan(twist.dot(n_k.cross(t_k)) / (twist.dot(n_k)));
    }

}   // namespace DCurvenet