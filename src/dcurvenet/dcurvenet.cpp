#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <cmath>
#include <utility>

namespace DCurvenet {
    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(Curvenet::curvenet* CN, int alpha): CN(CN) {
        std::vector<Curvenet::Control> cnCtrl = CN->controls();
        std::vector<Curvenet::HalfEdge> cnHE = CN->halfedges();
        std::vector<Curvenet::CubicSpline> cnSpline = CN->splines();
        std::vector<Curvenet::Curve> cnCrv = CN->curves();
        // Defensive reset
        V.clear();
        E.clear();
        HE.clear();
        C.clear();
        // 1. First copy in the control points
        int num_controls = cnCtrl.size();
        for (int c = 0; c < num_controls; c++) {
            int v = V.size();
            V.emplace_back();
            V[v].n = cnCtrl[c].n;
            V[v].pos = cnCtrl[c].pos;
            V[v].cn_idx = c;
            V[v].cn_type = cnCtrl[c].cType;
        }
        // 2. Initialize new dCN copies of the CN curves
        for (int crv = 0; crv < cnCrv.size(); crv++) {
            const std::vector<int>& crvSplines = cnCrv[crv].splines;
            int c = C.size();
            C.emplace_back();
            C[c].cn_idx = crv;
            // Track the halfedges between splines in the curve
            int prev_he0 = -1;
            int next_he1 = -1;
            int temp_origin = -1;   // Tracks the new halfedge's origin
            // For each curve, initialize all its splines (+ edges, halfedges, verts)
            for (int s_idx = 0; s_idx < crvSplines.size(); s_idx++) {
                int s = crvSplines[s_idx];
                // Get start and endpoint controls
                int s_start = cnHE[cnSpline[s].he].origin;
                int s_end = cnHE[cnHE[cnSpline[s].he].twin].origin;
                if (s_idx == 0) {   // Initialize curve start on first spline
                    C[c].start = s_start;
                    temp_origin = s_start;
                } else if (s_idx == crvSplines.size() - 1) {
                    C[c].end = s_end;
                }

                // TODO: Compute actual number of samples to use
                // Compute samples
                double arclen = CN->arclenEst(s);
                int n_samples = std::max(3, 50);
                // Assume n_samples will always be > 3
                std::vector<Eigen::Vector3d> samples = CN->unifSample(s, n_samples);
                // Get local spline idx. If multiple, then this spline must make a self-loop
                std::vector<int> startLocalSplineIdx = CN->controlLocalSplineIdx(s, s_start);
                std::vector<int> endLocalSplineIdx = CN->controlLocalSplineIdx(s, s_end);
                if (startLocalSplineIdx.size() == 2) {  // self-loop
                    endLocalSplineIdx[0] = startLocalSplineIdx[1];
                }
                // Start outgoing edge
                for (int i = 1; i < n_samples; i++) {
                    int v;
                    // Initialize new structures
                    if (i == n_samples -1) {    // Last sample
                        v = s_end;
                    } else {
                        v = V.size();
                        V.emplace_back();
                    }
                    int he0 = HE.size();
                    int he1 = he0+1;
                    HE.emplace_back();
                    HE.emplace_back();
                    int e = E.size();
                    E.emplace_back();
                    
                    // Rewire
                    Eigen::Vector3d edgeVec = V[v].pos - V[temp_origin].pos;
                    HE[he0].twin = he1;
                    HE[he1].twin = he0;
                    HE[he0].prev = prev_he0;
                    HE[he1].next = next_he1;
                    HE[he0].tangent = edgeVec.normalized();
                    HE[he1].tangent = -1 * edgeVec.normalized();
                    HE[he0].l = edgeVec.norm();
                    HE[he1].l = edgeVec.norm();
                    if (i == n_samples - 1) {    // End vertex is the next
                        V[v].adjHE[endLocalSplineIdx[0]] = he1;  // Add he1 to outgoing of end
                        if (s_idx == crvSplines.size()-1) {
                            C[c].he_end = he1;
                        }
                    } else if (i == 1) { // Start vertex is the prev
                        V[v].adjHE[startLocalSplineIdx[0]] = he0;  // Add he0 to outgoing
                        // Also set the first halfedge of the curve
                        if (s_idx == 0) {
                            C[c].he_start = he0;
                        }
                    } else { // Middle vertex (symmetrize next and prev)
                        HE[prev_he0].next = he0;
                        HE[next_he1].prev = he1;
                        V[v].adjHE.push_back(he0);
                    }
                    HE[he0].dest = v;
                    HE[he1].dest = temp_origin;
                    HE[he0].edge = e;
                    HE[he1].edge = e;
                    HE[he0].sign = true;    // left side
                    HE[he1].sign = false;   // right side
                    E[e].he = he0;
                    E[e].curve = c;
                    // Update loop params
                    prev_he0 = he0;
                    next_he1 = he1;
                    temp_origin = v;
                }
            }
        }

        // 3. Iterate over every control vertex that is an intersection/anchor and rewire its outgoing/incoming halfedges
        for (int v = 0; v < num_controls; v++) {
            if (V[v].cn_type == 2) {
                continue;
            }
            const std::vector<int> adjHE = V[v].adjHE;
            for (int he = 0; he < adjHE.size(); he++) {
                int he0 = adjHE[he];
                int he1 = adjHE[(he0+1)%adjHE.size()];
                HE[he0].prev = HE[he1].twin;
                HE[HE[he1].twin].next = he0;
            }
        }

        // 4. Compute all corner normals and widths
        allCornerNormalsAndWidths();
        // 5. Transport normals and widths along all splines
        transportNormalsAndWidths();
        // 6. Compute scaled frames on all splines
        computeScaledFrames();
    }

    // For intersections, computes their corner normals. For non-controls, this method does nothing (return -1)
    int dcurvenet::vertCornerNormalsWidths(int v) {
        // Do not process vertices which are not marked as intersections or anchors
        if ((V[v].cn_type != 1 && V[v].cn_type != 3) || !V[v].active) {
            return -1;
        }
        double eps = 1e-6;
        const std::vector<int> adjHE = V[v].adjHE;
        std::vector<Eigen::Vector3d> adjNormals(adjHE.size());
        // Explictly handle anchors -- but only if the outgoing curve has an endpoint anchor
        if (V[v].cn_type == 1) {
            // Grab the adjacent curve
            int he = adjHE[0];
            int c = E[HE[he].edge].curve;
            // Figure out if v is at the end or start of the curve + assign normals and widths
            if (HE[he].sign) {
                C[c].N_pos.first = V[v].n;
                C[c].W_pos.first = HE[he].l;
                C[c].N_neg.first = V[v].n;
                C[c].W_neg.first = HE[he].l;
            } else {
                C[c].N_pos.second = V[v].n;
                C[c].W_pos.second = HE[he].l;
                C[c].N_neg.second = V[v].n;
                C[c].W_neg.second = HE[he].l;
            }
            // Note that if the other endpoint of the curve is an intersection, these will be ignored.
            return 1;
        }
        // Otherwise, we need to compute normals explicitly
        std::vector<bool> skipList(adjNormals.size(), false);
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he0+1)%adjHE.size()];
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;
            // First check if we are parallel. If so skip for now
            double dotProdTest = HE[he0].tangent.dot(HE[he1].tangent);
            if (dotProdTest >= 1.0-eps || dotProdTest <= -1.0+eps) {
                skipList[he] = true;
                continue;
            }
            
            // Compute corner normal
            Eigen::Vector3d cornerNormal = HE[he0].tangent.cross(HE[he1].tangent);

            // Assign corner normals to curves
            // First figure out if this is the start halfedge of the curve
            if (HE[he0].sign) {
                C[c0].N_pos.first = cornerNormal.normalized();
            } else {
                C[c0].N_neg.second = cornerNormal.normalized();
            }
            if (HE[he1].sign) {
                C[c1].N_neg.first = cornerNormal.normalized();
            } else {
                C[c1].N_pos.second = cornerNormal.normalized();
            }
            adjNormals[he] = cornerNormal;
        }
        // Loop over any degenerate cases (straight angles)
        for (int he = 0; he < skipList.size(); he++) {
            int he0 = adjHE[he];
            int he_m1 = adjHE[(he + adjHE.size() - 1)%adjHE.size()];
            int he1 = adjHE[(he + adjHE.size() - 1)%adjHE.size()];
            // Get adjacent curves
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;

            // Average the adjacent vectors
            Eigen::Vector3d cornerNormal = (adjNormals[he_m1] + adjNormals[he1])/((adjNormals[he_m1] + adjNormals[he1]).norm());
            // Assign corner normals to curves
            if (HE[he0].sign) {
                C[c0].N_pos.first = cornerNormal;
            } else {
                C[c0].N_neg.second = cornerNormal;
            }
            if (HE[he1].sign) {
                C[c1].N_neg.first = cornerNormal;
            } else {
                C[c1].N_pos.second = cornerNormal;
            }
            adjNormals[he] = cornerNormal;
        }
        // Compute widths
        // Use the computed normal norm list to calculate widths
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he0+1)%adjHE.size()];
            int he_m1 = adjHE[(he0-1+adjHE.size())%adjHE.size()];
            int c0 = E[HE[he0].edge].curve;
            double cornerNormalNorm0;
            double he0_len = HE[he0].l;
            double he_m1_len = HE[he1].l;
            double he1_len = HE[he1].l;

            // Compute corner widths
            double cornerWidth0 = he0_len + adjNormals[he].norm()*(he1_len - he0_len);
            double cornerWidth1 = he0_len + adjNormals[he_m1].norm()*(he_m1_len - he0_len);

            // Assign corner normals to curves
            // First figure out if this is the start halfedge of the curve
            if (HE[he0].sign) {
                C[c0].W_pos.first = cornerWidth0;
                C[c0].W_neg.first = cornerWidth1;
            } else {
                C[c0].W_neg.second = cornerWidth0;
                C[c0].W_pos.second = cornerWidth1;
            }
        }
        return 1;
    }

    // Corner normals on only intersections
    int dcurvenet::allCornerNormalsAndWidths() {
        for (int v = 0; v < V.size(); v++) {
            vertCornerNormalsWidths(v);
        }
        return 1;
    }

    // Transport corner normals and widths from the two ends of a curve
    int dcurvenet::transportNWOnCurve(int c) {
        if (!C[c].active) {
            return -1;
        }
        // First, query what kinds of endpoints we have
        int start = C[c].start;
        int end = C[c].end;
        int he_start = C[c].he_start;
        int he_end = C[c].he_end;
        int start_type = V[C[c].start].cn_type;
        int end_type = V[C[c].end].cn_type;

        // We haven't set the normals and widths for closed curves yet
        if (start_type == 2 && end_type == 2) {
            C[c].N_pos.first = V[start].n;
            C[c].N_pos.second = V[start].n;
            C[c].N_neg.first = V[start].n;
            C[c].N_neg.second = V[start].n;

            C[c].W_pos.first = HE[he_start].l;
            C[c].W_pos.second = HE[he_end].l;
            C[c].W_neg.first = HE[he_start].l;
            C[c].W_neg.second = HE[he_end].l;
        }

        // Case 1: Both endpoints are intersections or both are anchors
        if (start_type == end_type) {
            // POSITIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = C[c].N_pos.first;
            Eigen::Vector3d end_n = C[c].N_pos.second;
            double start_w = C[c].W_pos.first;
            double end_w = C[c].W_pos.second;
            // 1. For the pos side, first trace until the end vertex, accumulating rotations and alpha values
            double total_len = accumulateRotations(C[c].he_start, end, rots, lens);
            // 2. Compute the torsion angle theta
            double torsion = computeTorsion(start_n, end_n, rots[rots.size()-1], -1*HE[C[c].he_end].tangent);
            // 3. Propagate rotations and widths to halfedges
            int he_curr = he_start;
            for (int he = 0; he < rots.size(); he++) {
                double alpha = lens[he]/total_len;
                Eigen::Matrix3d he_torsion = Utils::computeRotation(HE[he_curr].tangent, alpha*torsion);
                HE[he_curr].normal = he_torsion * rots[he] * start_n;
                HE[he_curr].w = (1-alpha)*start_w + (alpha)*end_w;
                he_curr = HE[he_curr].next;
            }

            // NEGATIVE SIDE
            // The same thing but backwards (so we can trace using the same function of next halfedges)
            start_n = C[c].N_neg.second;
            end_n = C[c].N_neg.first;
            start_w = C[c].W_neg.second;
            end_w = C[c].W_neg.first;
            total_len = accumulateRotations(C[c].he_end, start, rots, lens);
            torsion = computeTorsion(start_n, end_n, rots[rots.size()-1], -1*HE[C[c].he_start].tangent);
            he_curr = he_end;
            for (int he = 0; he < rots.size(); he++) {
                double alpha = lens[he]/total_len;
                Eigen::Matrix3d he_torsion = Utils::computeRotation(HE[he_curr].tangent, alpha*torsion);
                HE[he_curr].normal = he_torsion * rots[he] * start_n;
                HE[he_curr].w = (1-alpha)*start_w + (alpha)*end_w;
                he_curr = HE[he_curr].next;
            }
            return 1;
        }
        // Case 2: Start is an intersection, end is an anchor
        else if ((start_type == 3) && (end_type == 1)) {
            // POSITIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = C[c].N_pos.first;
            double start_w = C[c].W_pos.first;
            // 1. For the pos side, first trace until the end vertex, accumulating rotations and alpha values
            double total_len = accumulateRotations(C[c].he_start, end, rots, lens);
            // 3. Propagate rotations and widths to halfedges
            int he_curr = he_start;
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].normal = rots[he] * start_n;
                HE[he_curr].w = start_w;
                he_curr = HE[he_curr].next;
            }

            // NEGATIVE SIDE
            start_n = C[c].N_neg.first;
            start_w = C[c].W_neg.first;
            he_curr = HE[he_start].twin;
            // Do not re-initialize rotations, since we only trace from interesection to anchor
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].normal = rots[he] * start_n;
                HE[he_curr].w = start_w;
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
            Eigen::Vector3d start_n = C[c].N_neg.second;
            double start_w = C[c].W_neg.second;
            int total_len = accumulateRotations(C[c].he_end, start, rots, lens);
            int he_curr = he_end;
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].normal = rots[he] * start_n;
                HE[he_curr].w = start_w;
                he_curr = HE[he_curr].next;
            }

            start_n = C[c].N_pos.second;
            start_w = C[c].W_pos.second;
            he_curr = HE[he_end].twin;
            // Do not re-initialize rotations, since we only trace from interesection to anchor
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].normal = rots[he] * start_n;
                HE[he_curr].w = start_w;
                he_curr = HE[he_curr].prev;
            }
            return 1;
        }
        return -1;
    }
    // Transport normals for all curves
    int dcurvenet::transportNormalsAndWidths() {
        for (int c = 0; c < C.size(); c++) {
            transportNWOnCurve(c);
        }
        return 1;
    }

    // Complete the scaled frames on a curve by computing the binormal and height
    int dcurvenet::computeScaledFrameOnCurve(int c) {
        if (!C[c].active) {
            return -1;
        }
        int start = C[c].start;
        int end = C[c].end;
        int he_start = C[c].he_start;
        int he_end = C[c].he_end;
        
        int he_curr = -1;
        // Iterate over every halfedge in the curve
        do {
            if (he_curr == -1) {
                he_curr = he_start;
            } else {
                he_curr = HE[he_curr].next;
            }
            // Positive side
            HE[he_curr].binormal = HE[he_curr].tangent.cross(HE[he_curr].normal);
            HE[he_curr].h = std::sqrt(HE[he_curr].l * HE[he_curr].w);
            // Negative side
            int he_neg = HE[he_curr].twin;
            HE[he_neg].binormal = HE[he_neg].tangent.cross(HE[he_neg].normal);
            HE[he_neg].h = std::sqrt(HE[he_neg].l * HE[he_neg].w);
        } while (HE[he_curr].dest != end && he_curr != -1);
        return 1;
    }
    int dcurvenet::computeScaledFrames() {
        for (int c = 0; c < C.size(); c++) {
            computeScaledFrameOnCurve(c);
        }
        return 1;
    }

}   // namespace DCurvenet