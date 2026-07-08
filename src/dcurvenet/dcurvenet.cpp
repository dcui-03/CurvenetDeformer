#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <cmath>
#include <utility>

namespace DCurvenet {
    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(Curvenet::curvenet* CN, double meanE, int alpha, int mode): CN(CN), alpha(alpha), meanE(meanE) {
        const std::vector<Curvenet::Control>& cnCtrl = CN->controls();
        int num_curves = CN->numCurves();
        // Defensive reset
        V.clear();
        HE.clear();
        C.clear();
        // 1. First copy in the control points
        int num_controls = cnCtrl.size();
        std::vector<int> ctrlVerts(num_controls);
        for (int c = 0; c < num_controls; c++) {
            int new_v = addVert(cnCtrl[c], c);
            inputCtoV[c] = new_v;
            ctrlVerts[c] = new_v;
        }
        // 2. Initialize new dCN copies of the CN curves
        for (int crv = 0; crv < num_curves; crv++) {
            int c = addCurve(crv);
            inputCrvToC[crv] = c;
        }

        // 3. Iterate over every control vertex that is an intersection/anchor and rewire its outgoing/incoming halfedges
        // NOTE: This is not really necessary. We already have the outgoing order of halfedges per control vertex, and we will never need to traverse betweem curves
        for (int v_idx = 0; v_idx < ctrlVerts.size(); v_idx++) {
            int v = ctrlVerts[v_idx];
            rewireVertAdjHE(v);
        }

        if (mode == 0) {
            // 4. Compute all corner normals and widths
            allCornerNormalsAndWidths();
            // 5. Transport normals and widths along all splines
            transportNormalsAndWidths();
            // 6. Compute scaled frames on all splines
            computeScaledFrames();

            // 7. Cleanup by copying realtime frames to rest frames
            copyFramesToNew();
        }
    }

    // Update the new frames on all halfedges
    void dcurvenet::updateDiscCurveNet() {
        const std::vector<Curvenet::Control>& cnCtrl = CN->C;
        const std::vector<Curvenet::CubicSpline>& cnSpline = CN->S;
        const std::vector<Curvenet::Curve>& cnCurve = CN->Crv;
        // 1. Copy new control positions to their corresponding dvert
        for (const auto& idxPair : inputCtoV) {
            V[idxPair.second].new_pos = cnCtrl[idxPair.first].new_pos;
        }
        // 2. Iterate over curves and recompute
        for (const auto& idxPair : inputCrvToC) {
            std::vector<int> splines = cnCurve[idxPair.first].splines;
            int c = idxPair.second;
            int he_curr = C[c].he_start;
            int v_prev = C[c].start;
            for (int s_idx = 0; s_idx < splines.size(); s_idx++) {
                int s = splines[s_idx];
                int n_samples = cnSpline[s].num_samples;
                // Assume n_samples will always be > 3
                const std::vector<Eigen::Vector3d> samples = CN->unifSample(s, n_samples);
                // Exploit the fact that this matches the halfedge direction that the curve was constructed from
                for (int i = 1; i < samples.size(); i++) {
                    // for numerical reasons, only copy in non-controls
                    if (i != samples.size() - 1) {
                        V[HE[he_curr].dest].new_pos = samples[i];
                    }
                    Eigen::Vector3d tangent = V[HE[he_curr].dest].new_pos - V[v_prev].new_pos;
                    HE[he_curr].tangent = tangent.normalized();
                    HE[he_curr].l = tangent.norm();
                    HE[HE[he_curr].twin].tangent = -1 * tangent.normalized();
                    HE[HE[he_curr].twin].l = tangent.norm();
                    he_curr = HE[he_curr].next;
                    v_prev = HE[he_curr].dest;
                }
            }
        }
        // Compute corner normals and widths
        allCornerNormalsAndWidths();
        // Transport all normals and widths
        transportNormalsAndWidths();
        // Compute scaled frames on all splines
        computeScaledFrames();
        return;
    }

    // Add a new vertex that matches an existing control
    int dcurvenet::addVert(Curvenet::Control ctrl, int ctrl_idx) {
        int v = V.size();
        V.emplace_back();
        V[v].n = ctrl.n;
        V[v].pos = ctrl.pos;
        V[v].new_pos = ctrl.pos;
        V[v].cn_idx = ctrl_idx;
        V[v].cn_type = ctrl.cType;
        V[v].adjHE.resize(ctrl.adjHE.size(), -1);
        return v;
    }
    // Add a vertex given its parameters
    int dcurvenet::addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n, int ctrl_idx, int ctrl_type, int adjSize) {
        int v = V.size();
        V.emplace_back();
        V[v].n = new_n;
        V[v].pos = new_pos;
        V[v].new_pos = new_pos;
        V[v].cn_idx = ctrl_idx;
        V[v].cn_type = ctrl_type;
        V[v].adjHE.resize(adjSize, -1);
        return v;
    }
    // Add a new edge in and return its halfedges
    int dcurvenet::addEdge(int origin, int dest, int prev_he0, int next_he1, int c) {
        if (origin > V.size() || dest > V.size()) {
            return -1;
        }
        int e = E.size();
        E.emplace_back();
        int he0 = HE.size();
        int he1 = he0+1;
        HE.emplace_back();
        HE.emplace_back();
        
        // Rewire
        E[e].he = he0;
        E[e].curve = c;
        HE[he0].dest = dest;
        HE[he1].dest = origin;
        Eigen::Vector3d edgeVec = V[dest].pos - V[origin].pos;
        HE[he0].twin = he1;
        HE[he1].twin = he0;
        HE[he0].prev = prev_he0;
        HE[he1].next = next_he1;
        if (prev_he0 != -1) {
            HE[prev_he0].next = he0;
        }
        if (next_he1 != -1) {
            HE[next_he1].prev = he1;
        }
        HE[he0].tangent = edgeVec.normalized();
        HE[he1].tangent = -1 * edgeVec.normalized();
        HE[he0].l = edgeVec.norm();
        HE[he1].l = edgeVec.norm();
        HE[he0].sign = true;    // left side
        HE[he1].sign = false;   // right side

        return e;
    }
    // Add a curve that matches an input curvenet curve
    int dcurvenet::addCurve(int crv) {
        const std::vector<Curvenet::HalfEdge>& cnHE = CN->halfedges();
        const std::vector<Curvenet::CubicSpline>& cnSpline = CN->splines();
        const std::vector<Curvenet::Curve>& cnCrv = CN->curves();

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
            } 
            if (s_idx == crvSplines.size() - 1) {   // Initialize curve end on the last spline
                C[c].end = s_end;
            }

            // Compute samples
            int n_samples = cnSpline[s].num_samples;
            // Assume n_samples will always be > 3
            const std::vector<Eigen::Vector3d> samples = CN->unifSample(s, n_samples);
            // Get local spline idx. If multiple, then this spline must make a self-loop
            std::vector<int> startLocalSplineIdx = CN->controlLocalSplineIdx(s_start, s);
            std::vector<int> endLocalSplineIdx = CN->controlLocalSplineIdx(s_end, s);
            if (s_start == s_end && startLocalSplineIdx.size() == 2) {  // Check self-loops
                endLocalSplineIdx[0] = startLocalSplineIdx[1];
            }
            // Start outgoing edge
            for (int i = 1; i < n_samples; i++) {
                int v;
                // Initialize new structures
                if (i == n_samples -1) {    // Last sample
                    v = s_end;
                } else {    // Initialize a new sample
                    v = addVert(samples[i]);
                }
                int new_edge = addEdge(temp_origin, v, prev_he0, next_he1, c);
                int he0 = E[new_edge].he;
                int he1 = HE[he0].twin;
                
                if (i == n_samples - 1) {    // End vertex is the next
                    V[v].adjHE[endLocalSplineIdx[0]] = he1;  // Add he1 to outgoing of end
                    if (s_idx == crvSplines.size()-1) {
                        C[c].he_end = he1;
                    }
                } else if (i == 1) { // Start vertex is the prev
                    V[temp_origin].adjHE[startLocalSplineIdx[0]] = he0;  // Add he0 to outgoing
                    // Also set the first halfedge of the curve
                    if (s_idx == 0) {
                        C[c].he_start = he0;
                    }
                } else { // Middle vertex (symmetrize next and prev)
                    V[temp_origin].adjHE.push_back(he0);
                }
                // Update loop params
                prev_he0 = he0;
                next_he1 = he1;
                temp_origin = v;
            }
        }
        return c;
    }
    // Rewire incoming and outgoing halfedges of intersection vertices
    int dcurvenet::rewireVertAdjHE(int v) {
        if (V[v].cn_type < 3) {
            return -1;
        }
        const std::vector<int> adjHE = V[v].adjHE;
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            HE[he0].prev = HE[he1].twin;
            HE[HE[he1].twin].next = he0;
        }
        return 1;
    }


    // For intersections, computes their corner normals. For non-controls, this method does nothing (return -1)
    int dcurvenet::vertCornerNormalsWidths(int v) {
        if (!V[v].active) {
            return -1;
        }
        double eps = 1e-6;
        const std::vector<int> adjHE = V[v].adjHE;
        std::vector<Eigen::Vector3d> adjNormals(adjHE.size());
        // Explictly handle anchors and closed curves
        if (V[v].cn_type < 3) {
            // Check if we are on the start of a curve
            int he0 = adjHE[0];
            int c = E[HE[he0].edge].curve;
            if ((C[c].start != v) && V[v].cn_type == 2) {   // valence 2 control that is not the start of the curve
                return -1;
            }
            for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
                int he = adjHE[he_idx];
                const Eigen::Vector3d& tan = HE[he].tangent;
                // Gram-schmidt to ensure the normal is orthogonal to each tangent
                Eigen::Vector3d n = V[v].n - V[v].n.dot(tan) * tan;
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
                double w = HE[he].l + l * (HE[next_he].l - HE[he].l);
                if (HE[he].sign) {
                    C[c].N_pos.first = n;
                    C[c].W_pos.first = w;
                    C[c].N_neg.first = n;
                    C[c].W_neg.first = w;
                } else {
                    C[c].N_pos.second = n;
                    C[c].W_pos.second = w;
                    C[c].N_neg.second = n;
                    C[c].W_neg.second = w;
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
            double dotProdTest = HE[he0].tangent.dot(HE[he1].tangent);
            if (dotProdTest <= -1.0+eps) {
                skipList[he] = true;
                continue;
            } else if (dotProdTest >= 1 - eps) {    // Two outgoing HEs are coincident
                const Eigen::Vector3d& tan = HE[he0].tangent;
                cornerNormal = V[v].n - V[v].n.dot(tan) * tan;
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
                cornerNormal = HE[he0].tangent.cross(HE[he1].tangent);
            }

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
            if (!skipList[he]) {
                continue;
            }
            int he0 = adjHE[he];
            int he_m1_local = (he + adjHE.size() - 1)%adjHE.size();
            int he_m1 = adjHE[he_m1];
            int he1_local = (he + 1)%adjHE.size();
            int he1 = adjHE[he1_local];
            // Get adjacent curves
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;

            // Average the adjacent vectors
            Eigen::Vector3d cornerNormal = adjNormals[he_m1_local] + adjNormals[he1_local];
            // Extremely unlikely, but just in case, put in a safeguard...
            if (cornerNormal.norm() <= eps) {
                Eigen::Vector3d tan = HE[he0].tangent;
                cornerNormal = V[v].n - V[v].n.dot(tan) * tan;cornerNormal = V[v].n - V[v].n.dot(tan) * tan;
                if (cornerNormal.norm() <= eps) {
                    // Just pick a random direction orthogonal to the tangent 
                    if (std::abs(tan(0)) < 0.9) {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
            }
            cornerNormal.normalize();

            // Assign corner normals to curves
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
        // Compute widths
        // Use the computed normal norm list to calculate widths
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            int he_m1 = adjHE[(he-1+adjHE.size())%adjHE.size()];
            int c0 = E[HE[he0].edge].curve;
            double cornerNormalNorm0;
            double he0_len = HE[he0].l;
            double he_m1_len = HE[he_m1].l;
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
            double total_len = accumulateRotations(C[c].he_end, start, rots, lens);
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
            // Re-orthogonalize normals for safety
            HE[he_curr].normal = (HE[he_curr].normal - HE[he_curr].normal.dot(HE[he_curr].tangent) * HE[he_curr].tangent).normalized();
            // Positive side
            HE[he_curr].binormal = (HE[he_curr].tangent.cross(HE[he_curr].normal)).normalized();
            HE[he_curr].h = std::sqrt(HE[he_curr].l * HE[he_curr].w);
            // Negative side
            int he_neg = HE[he_curr].twin;
            // Re-orthogonalize normals for safety
            HE[he_neg].normal = (HE[he_neg].normal - HE[he_neg].normal.dot(HE[he_neg].tangent) * HE[he_neg].tangent).normalized();
            HE[he_neg].binormal = (HE[he_neg].tangent.cross(HE[he_neg].normal)).normalized();
            HE[he_neg].h = std::sqrt(HE[he_neg].l * HE[he_neg].w);
        } while (he_curr != -1 && HE[he_curr].dest != end);
        return 1;
    }
    int dcurvenet::computeScaledFrames() {
        for (int c = 0; c < C.size(); c++) {
            computeScaledFrameOnCurve(c);
        }
        return 1;
    }

    // Copy scaled frame data to new local variables
    int dcurvenet::copyFrameToNew(int he) {
        if (he >= HE.size()) {
            return -1;
        }
        HE[he].rest_tangent = HE[he].tangent;
        HE[he].rest_binormal = HE[he].binormal;
        HE[he].rest_normal = HE[he].normal;

        HE[he].rest_l = HE[he].l;
        HE[he].rest_w = HE[he].w;
        HE[he].rest_h = HE[he].h;
        return 1;
    }
    int dcurvenet::copyFramesToNew() {
        for (int he = 0; he < HE.size(); he++) {
            copyFrameToNew(he);
        }
        return 1;
    }
}   // namespace DCurvenet