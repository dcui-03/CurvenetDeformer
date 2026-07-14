#include "cutmesh.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <utility>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <cmath>
#include <iostream>

namespace Mesh {

int cutmesh::polyscopeFormat(Eigen::MatrixXd& Verts, 
                            std::vector<std::vector<int>>& Faces) {
    if (active_v < 3) {
        return -1;
    }
    Verts.resize(active_v, 3);
    Faces.resize(active_f);
    std::map<int, int> vToPSV;
    int curr_v = 0;
    for (int v = 0; v < V.size(); v++) {
        if (V[v].active) {
            if (curr_v >= active_v) {
                std::cout << "Counted faces does not equal active faces" << std::endl;
                return -1;
            }
            Verts.row(curr_v) = V[v].pos;
            vToPSV[v] = curr_v;
            curr_v++;
        }
    }
    int curr_f = 0;
    for (int f = 0; f < F.size(); f++) {
        if (F[f].active) {
            if (curr_f >= active_f) {
                std::cout << "Counted faces does not equal active faces" << std::endl;
                return -1;
            }
            std::vector<int> adjV = faceAdjVertIdxs(f);
            for (int v = 0; v < adjV.size(); v++) {
                adjV[v] = vToPSV[adjV[v]];
            }
            Faces[curr_f] = adjV;
            curr_f++;
        }
    }
    return 1;
}


}   // namespace Mesh