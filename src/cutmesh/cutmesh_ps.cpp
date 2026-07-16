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
                            std::vector<std::vector<int>>& Faces, 
                            std::vector<glm::vec3>& VertN, 
                            std::vector<glm::vec3>& FaceN,
                            std::vector<glm::vec3>& cornerIdx,
                            std::vector<glm::vec3>& projVecs) const {
    if (active_v < 3) {
        return -1;
    }
    Verts.resize(active_v, 3);
    Faces.resize(active_f);
    VertN.resize(active_v);
    FaceN.resize(active_f);
    cornerIdx.resize(active_v);
    projVecs.resize(active_v);
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
            VertN[curr_v] = Utils::eigenToGLM(V[v].n);
            if (V[v].label == 1 || V[v].label == 2) {
                int c_idx = V[v].corner_idx;
                if (c_idx < 0) {
                    std::cout << "Cut vertex has a non-existent corner index" << std::endl;
                }
                Eigen::Vector3d edgeMidpoint = 0.5 * (V[HE[c_idx].dest].pos + V[HE[HE[c_idx].twin].dest].pos);
                cornerIdx[curr_v] = Utils::eigenToGLM(edgeMidpoint - V[v].pos);
                projVecs[curr_v] = Utils::eigenToGLM(V[v].defData.projVector);
            } else {
                cornerIdx[curr_v] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
                projVecs[curr_v] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
            }

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
                auto it = vToPSV.find(adjV[v]);
                if (it == vToPSV.end()) {
                    std::cout << "polyscopeFormat(): face references inactive/missing vertex "
                            << adjV[v] << " in face " << f << std::endl;
                    return -1;
                }
                adjV[v] = it->second;
            }
            Faces[curr_f] = adjV;
            FaceN[curr_f] = Utils::eigenToGLM(F[f].n);
            curr_f++;
        }
    }
    return 1;
}


}   // namespace Mesh