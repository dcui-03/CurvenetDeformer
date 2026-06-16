#include "curvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <cmath>
#include <array>


namespace Curvenet {
    // Initialize from an existing list of controls, splines
    curvenet::curvenet(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents, std::vector<std::array<int, 4>> Splines) {
        for (int c = 0; c < Controls.size(); c++) {
            addControl(Controls[c]);
        }
        for (int s = 0; s < Splines.size(); s++) {
            std::array<int, 4> S = Splines[s];
            addSpline(S[0], S[3], Tangents[S[1]], Tangents[S[2]]);
        }
        return;
    }

    // empty initializer
    curvenet::curvenet() {
        C.clear();
        HE.clear();
        S.clear();
    }

    // Create new control
    int curvenet::addControl(Eigen::Vector3d pos) {
        int c = C.size();
        C.emplace_back();
        C[c].pos = pos;
        return c;
    }
    // Add a spline to the spline list given indices of the points
    int curvenet::addSpline(int start, int end, Eigen::Vector3d t0, Eigen::Vector3d t1) {
        if (!C[start].active || !C[end].active) {
            return -1;
        }
        // Create 2 new halfedges and a new spline
        int he0 = HE.size();
        int he1 = he0+1;
        HE.emplace_back();
        HE.emplace_back();
        int s = S.size();
        S.emplace_back();
        // Re-wire
        HE[he0].twin = he1;
        HE[he1].twin = he0;
        HE[he0].origin = start;
        HE[he1].origin = end;
        HE[he0].tan = t0;
        HE[he1].tan = t1;
        S[s].he = he0;
        HE[he0].s = s;
        HE[he1].s = s;
        // Insert spline into vertex list
        C[start].adjHE.push_back(he0);
        C[end].adjHE.push_back(he1);
        // Re-sort the outgoing halfedges of each control
        //sortAdjHE(start);
        sortAdjHE(end);
        return s;
    }

    // Move control
    int curvenet::editControlPos(int c, const Eigen::Vector3d pos) {
        if (!C[c].active) {
            return -1;
        }
        C[c].pos = pos;
        return c;
    }
    // Change control normal
    int curvenet::editControlN(int c, Eigen::Vector3d normal) {
        // Catch degenerate cases
        if ((!C[c].active) || normal.squaredNorm() == 0.0) {
            return -1;
        }
        C[c].n = normal.normalized();   // Always normalize for safety
        return 1;
    }

    // Sort adjacent tangent vectors to a control point
    int curvenet::sortAdjHE(int c) {
        if ((!C[c].active) || (C[c].n == Eigen::Vector3d::Zero())) {
            return -1;
        }

        std::vector<int> adjHE = C[c].adjHE;
        std::vector<Eigen::Vector3d> adjT = ctrlAdjTans(c);

        // Compute a local angle for each adjacent halfedge
        std::vector<double> adjTheta;
        Eigen::Vector3d t1, t2;
        Utils::buildPlaneBasis(C[c].n, t1, t2);
        for (int t = 0; t < adjT.size(); t++) {
            Eigen::Vector3d target = C[c].pos + adjT[t];
            double theta;
            bool success = Utils::directionAngleInPlane(C[c].pos, target, C[c].n, t1, t2, theta);
            // Insert theta into the local adjacency list
            if (!success) { // If degenerate, just give it a big angle so that it gets sorted to the end
                theta = 3.0*M_PI;
            }
            adjTheta.push_back(theta);
        }

        // Sort by angle
        Utils::doubleListIdxSort(adjTheta, adjHE);
        // Do the necessary re-wiring
        for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
            int he = adjHE[he_idx];
            int he_next = adjHE[(he_idx + 1)%adjHE.size()];
            HE[he].prev = HE[he_next].twin;
            HE[HE[he_next].twin].next = he;
        }
        // Assign back to control
        C[c].adjHE = adjHE;
        C[c].sorted = true;
        return 1;   // success
    }
    int curvenet::sortAdjHEAll() {
        for (int c = 0; c < C.size(); c++) {
            if (C[c].active && !C[c].sorted) {
                sortAdjHE(c);
            }
        }
        return 1;
    }
    
    // Compute what kind of vertex each control is using the valence of splines
    int curvenet::assignCtrlType(int c) {
        std::vector<int> adjS = ctrlAdjSplines(c);
        return adjS.size();
    }
    int curvenet::assignCtrlTypeAll(int c) {
        for (int c = 0; c < C.size(); c++) {
            if (C[c].active) {
                assignCtrlType(c);
            }
        }
        return 1;
    }

}   // namespace Curvenet