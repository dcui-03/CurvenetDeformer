#include "profilemover.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include <stdexcept>

// Scaled frame / deformation gradient computation on the discrete curvenet (dCN).
// Moved here from dcurvenet, since this is ProfileMover-specific machinery.

namespace ProfileMover {

    void profilemover::initDeformation() {
        heDefData.resize(dCN.numHalfedges());
        computeEdgeFrames();
        std::vector<curveDeformData> curveData;
        allCornerNormalsAndWidths(curveData);
        transportNormalsAndWidths(curveData);
        computeScaledFrames();
        copyFramesToRest();
    }

    void profilemover::updateDeformation() {
        dCN.updateDiscCurveNet();
        computeEdgeFrames();
        std::vector<curveDeformData> curveData;
        allCornerNormalsAndWidths(curveData);
        transportNormalsAndWidths(curveData);
        computeScaledFrames();
        validateFrames();
        computeDefGradAll();
    }

    // Compute the raw edge tangent/length for every dCN halfedge from its current positions
    int profilemover::computeEdgeFrames() {
        std::vector<DCurvenet::Vert>& dCN_V = dCN.V;
        std::vector<DCurvenet::HalfEdge>& HE = dCN.HE;
        for (int he = 0; he < HE.size(); he++) {
            int origin = HE[HE[he].twin].dest;
            int dest = HE[he].dest;
            Eigen::Vector3d edgeVec = dCN_V[dest].new_pos - dCN_V[origin].new_pos;
            heDefData[he].newFrame.tangent = edgeVec.normalized();
            heDefData[he].newFrame.l = edgeVec.norm();
        }
        return 1;
    }

    // For intersections, computes their corner normals. For non-controls, this method does nothing (return -1)
    int profilemover::vertCornerNormalsWidths(int v, std::vector<curveDeformData>& curveData) {
        std::vector<DCurvenet::Vert>& dCN_V = dCN.V;
        std::vector<DCurvenet::HalfEdge>& HE = dCN.HE;
        std::vector<DCurvenet::Edge>& E = dCN.E;
        std::vector<DCurvenet::Curve>& dCN_C = dCN.C;
        if (!dCN_V[v].active) {
            return -1;
        }
        // Only original curvenet controls have corner data
        if (dCN_V[v].cn_idx < 0 || dCN_V[v].cn_type < 1) {
            return -1;
        }
        double eps = 1e-12;
        const std::vector<int> adjHE = dCN_V[v].adjHE;
        std::vector<Eigen::Vector3d> adjNormals(adjHE.size());
        // Explictly handle anchors and closed curves
        if (dCN_V[v].cn_type < 3) {
            // Check if we are on the start of a curve
            int he0 = adjHE[0];
            int c = E[HE[he0].edge].curve;
            if ((dCN_C[c].start != v) && dCN_V[v].cn_type == 2) {   // valence 2 control that is not the start of the curve
                return -1;
            }
            for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
                int he = adjHE[he_idx];
                const Eigen::Vector3d& tan = heDefData[he].newFrame.tangent;
                // Gram-schmidt to ensure the normal is orthogonal to each tangent
                Eigen::Vector3d n = dCN_V[v].n - dCN_V[v].n.dot(tan) * tan;
                if (n.norm() <= eps) {  // Extremely rare case where normal == tangent or normal == -tangent
                    // Just pick a random direction orthogonal to the tangent
                    if (std::abs(tan(0)) < 0.9) {
                        n = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        n = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
                double l = n.norm();
                n.normalize();
                int next_he = adjHE[(he_idx+1)%adjHE.size()];
                double w = heDefData[he].newFrame.l + l * (heDefData[next_he].newFrame.l - heDefData[he].newFrame.l);
                if (dCN.isPositiveHalfedge(he)) {
                    curveData[c].N_pos.first = n;
                    curveData[c].W_pos.first = w;
                    curveData[c].N_neg.first = n;
                    curveData[c].W_neg.first = w;
                } else {
                    curveData[c].N_pos.second = n;
                    curveData[c].W_pos.second = w;
                    curveData[c].N_neg.second = n;
                    curveData[c].W_neg.second = w;
                }
            }
            return 1;
        }

        // Otherwise, we are at an intersection and need to compute normals explicitly
        std::vector<bool> skipList(adjNormals.size(), false);
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;
            Eigen::Vector3d cornerNormal;
            // First check if we are parallel. If so skip for now
            double dotProdTest = (heDefData[he0].newFrame.tangent.normalized()).dot(heDefData[he1].newFrame.tangent.normalized());
            if (dotProdTest <= -1.0+eps) {
                skipList[he] = true;
                continue;
            } else if (dotProdTest >= 1 - eps) {    // Two outgoing HEs are coincident
                const Eigen::Vector3d& tan = heDefData[he0].newFrame.tangent;
                cornerNormal = dCN_V[v].n - dCN_V[v].n.dot(tan) * tan;
                if (cornerNormal.norm() <= eps) {  // Extremely rare case where normal == tangent or normal == -tangent
                    // Just pick a random direction orthogonal to the tangent
                    if (std::abs(tan(0)) < 0.9) {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
            } else {
                // Compute corner normal as usual
                cornerNormal = heDefData[he0].newFrame.tangent.cross(heDefData[he1].newFrame.tangent);
            }
            // Flip normals if the signed angle was > 180 (ex., simplified check against the vertex normal)
            if (cornerNormal.normalized().dot(dCN_V[v].n) < 0.0) {
                cornerNormal *= -1.0;
            }

            // Assign corner normals to curves
            // First figure out if this is the start halfedge of the curve
            if (dCN.isPositiveHalfedge(he0)) {
                curveData[c0].N_pos.first = cornerNormal.normalized();
            } else {
                curveData[c0].N_neg.second = cornerNormal.normalized();
            }
            if (dCN.isPositiveHalfedge(he1)) {
                curveData[c1].N_neg.first = cornerNormal.normalized();
            } else {
                curveData[c1].N_pos.second = cornerNormal.normalized();
            }
            adjNormals[he] = cornerNormal;
        }
        // Loop over any degenerate cases (straight angles)
        for (int he = 0; he < skipList.size(); he++) {
            if (!skipList[he]) {
                continue;
            }
            int he0 = adjHE[he];
            int he_m1_local = (he + adjHE.size() - 1)%adjHE.size();
            int he1_local = (he + 1)%adjHE.size();
            int he1 = adjHE[he1_local];
            // Get adjacent curves
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;

            // Average the adjacent vectors
            Eigen::Vector3d cornerNormal = adjNormals[he_m1_local] + adjNormals[he1_local];
            // Extremely unlikely, but just in case, put in a safeguard...
            if (cornerNormal.norm() <= eps) {
                Eigen::Vector3d tan = heDefData[he0].newFrame.tangent;
                cornerNormal = dCN_V[v].n - dCN_V[v].n.dot(tan) * tan;
                if (cornerNormal.norm() <= eps) {
                    // Just pick a random direction orthogonal to the tangent
                    if (std::abs(tan(0)) < 0.9) {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
            }

            // Assign corner normals to curves
            if (dCN.isPositiveHalfedge(he0)) {
                curveData[c0].N_pos.first = cornerNormal.normalized();
            } else {
                curveData[c0].N_neg.second = cornerNormal.normalized();
            }
            if (dCN.isPositiveHalfedge(he1)) {
                curveData[c1].N_neg.first = cornerNormal.normalized();
            } else {
                curveData[c1].N_pos.second = cornerNormal.normalized();
            }
            adjNormals[he] = cornerNormal;
        }
        // Compute widths
        // Use the computed normal norm list to calculate widths
        for (int he = 0; he < adjHE.size(); he++) {
            // Loop's halfedge indices
            int he0_local = he;
            int he1_local = (he + 1) % adjHE.size();
            int he_m1_local = (he - 1 + adjHE.size()) % adjHE.size();
            // Get the global halfedge index
            int he0 = adjHE[he0_local];
            int he1 = adjHE[he1_local];
            int he_m1 = adjHE[he_m1_local];
            // Get the curve and the corresponding lengths
            int c0 = E[HE[he0].edge].curve;
            double he0_len = heDefData[he0].newFrame.l;
            double he_m1_len = heDefData[he_m1].newFrame.l;
            double he1_len = heDefData[he1].newFrame.l;

            // Compute corner widths
            double cornerWidth0 = he0_len + adjNormals[he0_local].norm() * (he1_len - he0_len);
            double cornerWidth1 = he0_len + adjNormals[he_m1_local].norm() * (he_m1_len - he0_len);

            // Assign corner normals to curves
            // First figure out if this is the start halfedge of the curve
            if (dCN.isPositiveHalfedge(he0)) {
                curveData[c0].W_pos.first = cornerWidth0;
                curveData[c0].W_neg.first = cornerWidth1;
            } else {
                curveData[c0].W_neg.second = cornerWidth0;
                curveData[c0].W_pos.second = cornerWidth1;
            }
        }
        return 1;
    }

    // Corner normals on only intersections
    int profilemover::allCornerNormalsAndWidths(std::vector<curveDeformData>& curveData) {
        std::vector<DCurvenet::Vert>& dCN_V = dCN.V;
        std::vector<DCurvenet::Curve>& dCN_C = dCN.C;
        curveData.clear();
        curveData.resize(dCN_C.size());
        // TODO: Parallelize? NOTE: This may not be safe, since we operate on curveData simultaneously
        #pragma omp parallel for
        for (int v = 0; v < dCN_V.size(); v++) {
            vertCornerNormalsWidths(v, curveData);
        }
        return 1;
    }

    // Transport corner normals and widths from the two ends of a curve
    int profilemover::transportNWOnCurve(int c, const curveDeformData& cData) {
        std::vector<DCurvenet::Vert>& dCN_V = dCN.V;
        std::vector<DCurvenet::HalfEdge>& HE = dCN.HE;
        std::vector<DCurvenet::Curve>& dCN_C = dCN.C;
        if (!dCN_C[c].active) {
            return -1;
        }
        // First, query what kinds of endpoints we have
        int start = dCN_C[c].start;
        int end = dCN_C[c].end;
        int he_start = dCN_C[c].he_start;
        int he_end = dCN_C[c].he_end;
        int start_type = dCN_V[dCN_C[c].start].cn_type;
        int end_type = dCN_V[dCN_C[c].end].cn_type;

        // Case 1: Both endpoints are intersections or both are anchors
        if (start_type == end_type) {
            // POSITIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = cData.N_pos.first;
            Eigen::Vector3d end_n = cData.N_pos.second;
            double start_w = cData.W_pos.first;
            double end_w = cData.W_pos.second;
            // 1. For the pos side, first trace until the end vertex, accumulating rotations alpha values
            double total_len = accumulateRotations(dCN_C[c].he_start, end, rots, lens);
            // 2. Compute the torsion angle theta
            double torsion = computeTorsion(start_n, end_n, rots[rots.size()-1], -1*heDefData[dCN_C[c].he_end].newFrame.tangent);
            // 3. Propagate rotations and widths to halfedges
            int he_curr = he_start;
            for (int he = 0; he < rots.size(); he++) {
                double beta = lens[he]/total_len;
                Eigen::Matrix3d he_torsion = Utils::computeRotation(heDefData[he_curr].newFrame.tangent, beta*torsion);
                heDefData[he_curr].newFrame.normal = he_torsion * rots[he] * start_n;
                heDefData[he_curr].newFrame.w = (1-beta)*start_w + (beta)*end_w;
                he_curr = HE[he_curr].next;
            }

            // NEGATIVE SIDE
            // The same thing but backwards (so we can trace using the same function of next halfedges)
            start_n = cData.N_neg.second;
            end_n = cData.N_neg.first;
            start_w = cData.W_neg.second;
            end_w = cData.W_neg.first;
            total_len = accumulateRotations(dCN_C[c].he_end, start, rots, lens);
            torsion = computeTorsion(start_n, end_n, rots[rots.size()-1], -1*heDefData[dCN_C[c].he_start].newFrame.tangent);
            he_curr = he_end;
            for (int he = 0; he < rots.size(); he++) {
                double beta = lens[he]/total_len;
                Eigen::Matrix3d he_torsion = Utils::computeRotation(heDefData[he_curr].newFrame.tangent, beta*torsion);
                heDefData[he_curr].newFrame.normal = he_torsion * rots[he] * start_n;
                heDefData[he_curr].newFrame.w = (1-beta)*start_w + (beta)*end_w;
                he_curr = HE[he_curr].next;
            }
            return 1;
        }
        // Case 2: Start is an intersection, end is an anchor
        else if ((start_type == 3) && (end_type == 1)) {
            // POSITIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = cData.N_pos.first;
            double start_w = cData.W_pos.first;
            // 1. For the pos side, first trace until the end vertex, accumulating rotations and alpha values
            double total_len = accumulateRotations(dCN_C[c].he_start, end, rots, lens);
            // 3. Propagate rotations and widths to halfedges
            int he_curr = he_start;
            for (int he = 0; he < rots.size(); he++) {
                heDefData[he_curr].newFrame.normal = rots[he] * start_n;
                heDefData[he_curr].newFrame.w = start_w;
                he_curr = HE[he_curr].next;
            }

            // NEGATIVE SIDE
            start_n = cData.N_neg.first;
            start_w = cData.W_neg.first;
            he_curr = HE[he_start].twin;
            // Do not re-initialize rotations, since we only trace from interesection to anchor
            for (int he = 0; he < rots.size(); he++) {
                heDefData[he_curr].newFrame.normal = rots[he] * start_n;
                heDefData[he_curr].newFrame.w = start_w;
                he_curr = HE[he_curr].prev;
            }
            return 1;
        }
        // Case 3: Start is an anchor, end is an intersection
        else if ((start_type == 1) && (end_type == 3)) {
            // Same as case 2 but trace backwards instead of forwards
            // NEGATIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = cData.N_neg.second;
            double start_w = cData.W_neg.second;
            double total_len = accumulateRotations(dCN_C[c].he_end, start, rots, lens);
            int he_curr = he_end;
            for (int he = 0; he < rots.size(); he++) {
                heDefData[he_curr].newFrame.normal = rots[he] * start_n;
                heDefData[he_curr].newFrame.w = start_w;
                he_curr = HE[he_curr].next;
            }

            start_n = cData.N_pos.second;
            start_w = cData.W_pos.second;
            he_curr = HE[he_end].twin;
            // Do not re-initialize rotations, since we only trace from interesection to anchor
            for (int he = 0; he < rots.size(); he++) {
                heDefData[he_curr].newFrame.normal = rots[he] * start_n;
                heDefData[he_curr].newFrame.w = start_w;
                he_curr = HE[he_curr].prev;
            }
            return 1;
        }
        return -1;
    }
    // Transport normals for all curves
    int profilemover::transportNormalsAndWidths(const std::vector<curveDeformData>& curveData) {
        std::vector<DCurvenet::Curve>& dCN_C = dCN.C;
        // TODO: Parallelize?
        #pragma omp parallel for
        for (int c = 0; c < dCN_C.size(); c++) {
            transportNWOnCurve(c, curveData[c]);
        }
        return 1;
    }

    int profilemover::computeScaledFrames() {
        std::vector<DCurvenet::HalfEdge>& HE = dCN.HE;
        // TODO: Parallelize?
        #pragma omp parallel for
        for (int he = 0; he < HE.size(); he++) {
            if (HE[he].active) {
                // Re-orthogonalize normals for safety
                heDefData[he].newFrame.normal = (heDefData[he].newFrame.normal -
                                                    heDefData[he].newFrame.normal.dot(heDefData[he].newFrame.tangent) *
                                                    heDefData[he].newFrame.tangent).normalized();
                // Positive side
                heDefData[he].newFrame.binormal = (heDefData[he].newFrame.tangent.cross(heDefData[he].newFrame.normal)).normalized();
                heDefData[he].newFrame.h = std::sqrt(std::abs(heDefData[he].newFrame.l * heDefData[he].newFrame.w));
            }
        }
        return 1;
    }

    // Make sure frames don't
    int profilemover::validateFrames() {
        double eps = 1e-12;
        for (int he = 0; he < heDefData.size(); he++) {
            const scaledFrame& newFrame = heDefData[he].newFrame;
            if (newFrame.l <= eps || newFrame.w <= eps || newFrame.h <= eps) {
                throw std::runtime_error("validateFrames(): degenerate frame");
            }
            if (!std::isfinite(newFrame.l) || !std::isfinite(newFrame.w) || !std::isfinite(newFrame.h)) {
                throw std::runtime_error("validateFrames(): non-finite scaled frame scale");
            }
        }
        return 1;
    }
    // Copy scaled frame data to new local variables
    int profilemover::copyFrameToRest(int he) {
        if (he >= heDefData.size()) {
            return -1;
        }
        heDefData[he].restFrame.tangent = heDefData[he].newFrame.tangent;
        heDefData[he].restFrame.binormal = heDefData[he].newFrame.binormal;
        heDefData[he].restFrame.normal = heDefData[he].newFrame.normal;

        heDefData[he].restFrame.l = heDefData[he].newFrame.l;
        heDefData[he].restFrame.w = heDefData[he].newFrame.w;
        heDefData[he].restFrame.h = heDefData[he].newFrame.h;
        return 1;
    }
    int profilemover::copyFramesToRest() {
        for (int he = 0; he < heDefData.size(); he++) {
            copyFrameToRest(he);
        }
        return 1;
    }

    // Accumulate rotation matrices, starting from an initial halfedge and tracing forward until we hit the goal vertex
    double profilemover::accumulateRotations(int start_he, int end_v, std::vector<Eigen::Matrix3d>& rots, std::vector<double>& lens) {
        std::vector<DCurvenet::HalfEdge>& HE = dCN.HE;
        rots.clear();
        lens.clear();
        int he_curr = start_he;
        int he_prev = start_he;
        if (start_he < 0) {
            return 0.0;
        }
        Eigen::Vector3d tan_prev = heDefData[he_curr].newFrame.tangent;
        Eigen::Matrix3d curr_rot = Eigen::Matrix3d::Identity();
        double curr_len = 0.0;
        // Traverse halfedges until we hit the end vertex
        do {
            Eigen::Vector3d tan_curr = heDefData[he_curr].newFrame.tangent;
            curr_rot = Utils::computeRotation(tan_prev, tan_curr) * curr_rot;   // Left multiply to accumulate rotations
            rots.push_back(curr_rot);
            curr_len += heDefData[he_curr].newFrame.l;
            lens.push_back(curr_len);

            he_prev = he_curr;
            he_curr = HE[he_curr].next;
            tan_prev = tan_curr;
        } while ((HE[he_prev].dest != end_v) && (he_curr != -1));
        return curr_len;
    }

    // Total torsion
    double profilemover::computeTorsion(Eigen::Vector3d n_1, Eigen::Vector3d n_k, Eigen::Matrix3d Om_k, Eigen::Vector3d t_k) {
        Eigen::Vector3d twist = Om_k * n_1;
        double y = twist.dot(n_k.cross(t_k));
        double x = twist.dot(n_k);
        return std::atan2(y, x);
    }

    // Compute the deformation gradient of a halfedge
    Eigen::Matrix3d profilemover::computeHEDefGrad(int he) {
        const scaledFrame& heRestFrame = heDefData[he].restFrame;
        const scaledFrame& heNewFrame = heDefData[he].newFrame;
        Eigen::Matrix3d F = (heNewFrame.l / heRestFrame.l) * (heNewFrame.tangent * heRestFrame.tangent.transpose()) +
                            (heNewFrame.w / heRestFrame.w) * (heNewFrame.binormal * heRestFrame.binormal.transpose()) +
                            (heNewFrame.h / heRestFrame.h) * (heNewFrame.normal * heRestFrame.normal.transpose());
        return F;
    }

    int profilemover::computeDefGradAll() {
        for (int he = 0; he < heDefData.size(); he++) {
            heDefData[he].defGrad = computeHEDefGrad(he);
        }
        return 1;
    }

    // Pre-compute maps
    int profilemover::computedCNMats(Eigen::MatrixXd& f_dCN_flat, Eigen::MatrixXd& x_dCN) {
        std::vector<DCurvenet::Vert>& dCN_V = dCN.V;
        int num_HE = dCN.numHalfedges();
        int num_V  = dCN.numVerts();

        if (f_dCN_flat.rows() != num_HE || f_dCN_flat.cols() != 9) {
            f_dCN_flat.resize(num_HE, 9);
        }

        if (x_dCN.rows() != num_V || x_dCN.cols() != 3) {
            x_dCN.resize(num_V, 3);
        }

        #pragma omp parallel for
        for (int he = 0; he < num_HE; he++) {
            f_dCN_flat.row(he) =
                Utils::flattenMatrix3d(heDefData[he].defGrad).transpose();
        }

        #pragma omp parallel for
        for (int v = 0; v < num_V; v++) {
            x_dCN.row(v) = dCN_V[v].new_pos.transpose();
        }

        return 1;
    }

}   // namespace ProfileMover
