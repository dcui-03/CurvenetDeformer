#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <utility>

namespace DCurvenet {
    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(Curvenet::curvenet* CN, Mesh::mesh* M): CN(CN) {
        const std::vector<Curvenet::Control>& cnCtrl = CN->controls();
        int num_curves = CN->numCurves();
        // Defensive reset
        V.clear();
        HE.clear();
        E.clear();
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
        // Compute projection data
        computeProjData(M);
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
        // TODO: Can we parallelize this?
        // #pragma omp parallel for
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
                    int v_dest = HE[he_curr].dest;
                    v_prev = v_dest;
                    he_curr = HE[he_curr].next;
                }
            }
        }
        return;
    }

    dcurvenet::dcurvenet() {

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
        if (origin < 0 || origin >= V.size() || dest < 0 || dest >= V.size()) {
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
        HE[he0].twin = he1;
        HE[he1].twin = he0;
        HE[he0].prev = prev_he0;
        HE[he1].next = next_he1;
        HE[he0].edge = e;
        HE[he1].edge = e;
        if (prev_he0 != -1) {
            HE[prev_he0].next = he0;
        }
        if (next_he1 != -1) {
            HE[next_he1].prev = he1;
        }

        return e;
    }
    // Add a curve that matches an input curvenet curve
    int dcurvenet::addCurve(int crv) {
        const std::vector<Curvenet::HalfEdge>& cnHE = CN->HE;
        const std::vector<Curvenet::CubicSpline>& cnSpline = CN->S;
        const std::vector<Curvenet::Curve>& cnCrv = CN->Crv;
        if (crv < 0 || crv >= cnCrv.size()) {
            return -1;
        }
        const std::vector<int>& crvSplines = cnCrv[crv].splines;
        if (crvSplines.empty()) {
            return -1;
        }
        int c = C.size();
        C.emplace_back();
        C[c].cn_idx = crv;

        // Track the previous positive halfedge and previous negative halfedge
        // so that each newly added segment is wired into the curve chain.
        int prev_he0 = -1;
        int next_he1 = -1;
        // Current origin vertex in the discrete curvenet
        int temp_origin = -1;

        for (int s_idx = 0; s_idx < crvSplines.size(); s_idx++) {
            int s = crvSplines[s_idx];
            if (s < 0 || s >= cnSpline.size()) {
                return -1;
            }
            // Get the CN halfedges for this spline
            int cn_he = cnSpline[s].he;
            int cn_he_twin = cnHE[cn_he].twin;
            // Get the start and end CN control indices
            int cn_start = cnHE[cn_he].origin;
            int cn_end   = cnHE[cn_he_twin].origin;
            // Get the start and end dCN vertices
            int s_start = inputCtoV.at(cn_start);
            int s_end   = inputCtoV.at(cn_end);
            if (s_idx == 0) {
                C[c].start = s_start;
                temp_origin = s_start;
            }
            if (s_idx == static_cast<int>(crvSplines.size()) - 1) {
                C[c].end = s_end;
            }

            // Sample the spline
            int n_samples = cnSpline[s].num_samples;
            const std::vector<Eigen::Vector3d> samples = CN->unifSample(s, n_samples);
            // Get the local spline index for the start and end vertices
            std::vector<int> startLocalSplineIdx = CN->controlLocalSplineIdx(cn_start, s);
            std::vector<int> endLocalSplineIdx   = CN->controlLocalSplineIdx(cn_end, s);

            // Self-loop: parent curvenet has two local halfedges at the same control
            if (cn_start == cn_end && startLocalSplineIdx.size() == 2) {
                endLocalSplineIdx[0] = startLocalSplineIdx[1];
            }
            // Insert the discretized spline (i = 0 is the start vert)
            for (int i = 1; i < n_samples; i++) {
                int v = -1;

                if (i == n_samples - 1) {   // Last sample is the end vert
                    v = s_end;
                } else {    // Interior sample
                    v = addVert(samples[i]);
                }
                // Add a new edge
                int new_edge = addEdge(temp_origin, v, prev_he0, next_he1, c);

                int he0 = E[new_edge].he;       // positive/canonical direction
                int he1 = HE[he0].twin;         // negative/opposite direction

                // Assign adjacent halfedge(s) for each prev vert
                if (V[temp_origin].cn_idx < 0) {    // Interior vert
                    V[temp_origin].adjHE.clear();
                    V[temp_origin].adjHE.push_back(he0);
                } else if (i == 1) {    // Previous vert must be the spline start
                    V[temp_origin].adjHE[startLocalSplineIdx[0]] = he0;
                    if (s_idx == 0) {   // If we are starting the curve, make it the curve starting halfedge
                        C[c].he_start = he0;
                    }
                }
                // Assign outgoing halfedge for the end vert when the spline ends
                if (i == n_samples - 1) {
                    V[v].adjHE[endLocalSplineIdx[0]] = he1;
                    if (s_idx == crvSplines.size() - 1) {  // We are ending the curve, so add the outgoing halfedge to the end vert
                        C[c].he_end = he1;
                    }
                }

                prev_he0 = he0;
                next_he1 = he1;
                temp_origin = v;
            }
        }

        return c;
    }
    // Rewire incoming and outgoing halfedges of intersection and anchor vertices
    int dcurvenet::rewireVertAdjHE(int v) {
        const std::vector<int> adjHE = V[v].adjHE;
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            HE[he0].prev = HE[he1].twin;
            HE[HE[he1].twin].next = he0;
        }
        return 1;
    }


    // Compute projection data for this dcurvenet point
    int dcurvenet::computeProjData(Mesh::mesh* M) {
        for (int v = 0; v < V.size(); v++) {
            if (!V[v].active) {
                continue;
            }
            Mesh::meshBindData bindData;
            if (M->computeVBinding(V[v].pos, bindData) != 1) {
                return -1;
            }
            V[v].proj.coords = bindData.coords;
            V[v].proj.elType = bindData.elType;
            V[v].proj.elIdx = bindData.elIdx;
            V[v].proj.projVec = bindData.offset;
        }
        return 1;
    }

    // Propagate weights to rest of curvenet using Laplacian
    int dcurvenet::propagateWeights() {
        using T = Eigen::Triplet<double>;
        std::vector<T> tripletList;
        tripletList.reserve(E.size() * 4);  // Conservative overestimate
        // First, build Laplacian system
        Eigen::SparseMatrix<double> cnL(V.size(), V.size());
        // Eigen::VectorXd M(V.size());
        Eigen::VectorXd f(V.size());
        f.setZero();
        // M.setZero();
        std::vector<bool> fixed(V.size(), false);
        int num_fixed = 0;
        // Process fixed verts first
        for (int v = 0; v < V.size(); v++) {
            int v_type = V[v].cn_type;
            int v_idx = V[v].cn_idx;

            if (v_type != -1 && CN->C[v_idx].fixed_w) {
                fixed[v] = true;
                V[v].w = CN->C[v_idx].w;
                f[v] = V[v].w;
                tripletList.push_back(T(v, v, 1.0));
                num_fixed++;
            }
        }
        if (num_fixed == 0) {
            for (int v = 0; v < V.size(); v++) {
                V[v].w = 1.0;
            }
            return 1;
        }
        for (int e = 0; e < E.size(); e++) {
            int v0 = HE[E[e].he].dest;
            int v1 = HE[HE[E[e].he].twin].dest;

            int v0_idx = V[v0].cn_idx;
            int v1_idx = V[v1].cn_idx;
            // If both are fixed, then skip
            if (fixed[v0] && fixed[v1]) {   // Unlikely case, but check anyways
                continue;
            }
            double e_len = std::max((V[v1].pos - V[v0].pos).norm(), 1e-8);
            double weight = 1/e_len;
            bool f0 = fixed[v0];
            bool f1 = fixed[v1];
            // Add to diagonal entries
            if (!f0 && !f1) {
                tripletList.push_back(T(v0, v0,  weight));
                tripletList.push_back(T(v1, v1,  weight));
                tripletList.push_back(T(v0, v1, -weight));
                tripletList.push_back(T(v1, v0, -weight));
            } else if (!f0 && f1) {
                tripletList.push_back(T(v0, v0, weight));
                f[v0] += weight * V[v1].w;
            } else if (f0 && !f1) {
                tripletList.push_back(T(v1, v1, weight));
                f[v1] += weight * V[v0].w;
            }
            // M[v0] += e_len;
            // M[v1] += e_len;
        }
        cnL.setFromTriplets(tripletList.begin(), tripletList.end());
        // M *= 0.5;
        // Eigen::VectorXd RHS = M.asDiagonal() * f;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factorL;
        factorL.analyzePattern(cnL);
        factorL.factorize(cnL);
        Eigen::VectorXd new_weights = factorL.solve(f);
        
        // Redistribute weights
        for (int v = 0; v < V.size(); v++) {
            V[v].w = std::clamp(new_weights[v], 0.0, 1.0);
        }

        return 1;
    }
}   // namespace DCurvenet