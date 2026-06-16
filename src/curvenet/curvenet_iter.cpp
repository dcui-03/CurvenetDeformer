#include "curvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <algorithm>


namespace Curvenet {
    // Returns a CCW list of the tangents to a specified control
    std::vector<Eigen::Vector3d> curvenet::ctrlAdjTans(int c) {
        std::vector<Eigen::Vector3d> adjT;
        for (int he = 0; he < C[c].adjHE.size(); he++) {
            adjT.push_back(HE[C[c].adjHE[he]].tan);
        }
        return adjT;
    }

    // Returns a list of the splines adjacent to a control 
    // If no self loops, then these are CCW
    std::vector<int> curvenet::ctrlAdjSplines(int c) {
        std::vector<int> adjS;
        std::vector<int> adjHE = C[c].adjHE;
        for (int he = 0; he < adjHE.size(); he++) {
            if (std::find(adjS.begin(), adjS.end(), HE[adjHE[he]].s) == adjS.end()) {
                adjS.push_back(HE[adjHE[he]].s);
            }
        }
        return adjS;
    }

    // Returns vertices adjacent to a spline
    std::vector<int> curvenet::splineAdjCtrls(int s) {
        std::vector<int> adjC;
        int he0 = S[s].he;
        adjC.push_back(HE[he0].origin);   // First vertex
        if (HE[HE[he0].twin].origin != HE[he0].origin) {
            adjC.push_back(HE[HE[he0].twin].origin);
        }
        return adjC;
    }

}   // namespace Curvenet