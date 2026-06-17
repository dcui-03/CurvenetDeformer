#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include <Eigen/Core>
#include <vector>
#include <algorithm>

namespace DCurvenet {

    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(Curvenet::curvenet* CN, int alpha): CN(CN) {
        std::vector<Curvenet::Control> cnCtrl = CN->controls();
        std::vector<Curvenet::HalfEdge> cnHE = CN->halfedges();
        std::vector<Curvenet::CubicSpline> cnSpline = CN->splines();
        std::vector<Curvenet::Curve> cnCrv = CN->curves();
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
                // Insert into 
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
                    HE[he0].scale(0) = edgeVec.norm();
                    HE[he1].scale(0) = edgeVec.norm();
                    if (i == n_samples - 1) {    // End vertex is the next
                        V[v].adjHE[endLocalSplineIdx[0]] = he1;  // Add he0 to outgoing
                    } else if (i == 1) { // Start vertex is the prev
                        V[v].adjHE[startLocalSplineIdx[0]] = he0;  // Add he0 to outgoing
                        // Also set the first halfedge of the curve
                        C[c].he = he0;
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
        // Do not process vertices which are not marked as intersections
        if (V[v].cn_type < 3) {
            return -1;
        }
        double eps = 1e-6;
        const std::vector<int> adjHE = V[v].adjHE;
        std::vector<Eigen::Vector3d> adjNormals(adjHE.size());
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
            double he0_len = HE[he0].scale(0);
            double he_m1_len = HE[he1].scale(0);
            double he1_len = HE[he1].scale(0);

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
        // First, query what kinds of endpoints we have
        int start = C[c].start;
        int end = C[c].end;
        int start_type = V[C[c].start].cn_type;
        int end_type = V[C[c].end].cn_type;

        // TODO: Maybe write a function which does tracing and rotation/alpha accumulation for you, given a start vert index and end vert index
        // ex. pass in two vectors, one for alpha (std::vector<double>), one for rotation (std::vector<Eigen::Matrix3d>). returns total alpha as a double.
        // TODO: Need a utils function which converts an angle into a rotation matrix.

        // Case 1: Both endpoints are intersections
        if ((start_type == 3) && (end_type == 3)) {
            // 1. For the pos side, first trace until the end vertex, accumulating rotations and alpha values
            // 2. Compute the torsion angle theta
            // 3. Accumulate torsion and apply to curve

            // 3. Repeat for negative side (but starting from the end and tracing to the start)
        }
        // Case 2: Both endpoints are anchors
        else if ((start_type == 1) && (end_type == 1)) {
            // Need to handle this case somehow, so set all normals as 
            int he0 = adjHE[0];
            int c = E[HE[he0].edge].curve;
            C[c].N_pos.first = V[v].n;
            C[c].N_pos.second = V[v].n;
            C[c].N_neg.first = V[v].n;
            C[c].N_neg.second = V[v].n;
            return 1;
            // Take the start and end normals as their associated vertex normal
            // Take the widths as ???? just the length of the tangent?
            // Then, apply same logic as the intersection version
        }
        // Case 3: Start is an intersection, end is an anchor
        else if ((start_type == 3) && (end_type == 1)) {
            // 1. Accumulate rotations
            // 2. Propagate normals
        }
        // Case 4: Start is an anchor, end is an intersection
        else {
            // Same as case 3 but trace backwards instead of forwards
        }
        return 1;
    }
    // Transport normals for all curves
    int dcurvenet::transportNormalsAndWidths() {
        for (int c = 0; c < C.size(); c++) {
            transportNWOnCurve(c);
        }
        return 1;
    }

}   // namespace DCurvenet