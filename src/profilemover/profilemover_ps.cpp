#include "profilemover.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <map>
#include <glm/vec3.hpp>
#include <cmath>
#include <iostream>

// Polyscope viewer-formatting readers for dCN and CM. Moved here from
// dcurvenet/cutmesh since profilemover already keeps tabs on both, and the dCN
// reader needs the (now profilemover-owned) scaled frame data.

namespace ProfileMover {

int profilemover::dcurvenetPolyscopeFormat(Eigen::MatrixXd& Verts,
                            std::vector<std::array<int, 2>>& Edges,
                            std::vector<glm::vec3>& posEdgeTangents,
                            std::vector<glm::vec3>& posEdgeBinormals,
                            std::vector<glm::vec3>& posEdgeNormals,
                            std::vector<glm::vec3>& negEdgeTangents,
                            std::vector<glm::vec3>& negEdgeBinormals,
                            std::vector<glm::vec3>& negEdgeNormals,
                            std::vector<double>& weights) const {
    const std::vector<Polynet::Vert>& dCN_V = dCN.V;
    const std::vector<Polynet::HalfEdge>& HE = dCN.HE;
    const std::vector<Polynet::Edge>& E = dCN.E;
    Verts.resize(dCN_V.size(), 3);
    weights.resize(dCN_V.size());
    Edges.resize(E.size());
    posEdgeTangents.resize(E.size());
    posEdgeBinormals.resize(E.size());
    posEdgeNormals.resize(E.size());
    negEdgeTangents.resize(E.size());
    negEdgeBinormals.resize(E.size());
    negEdgeNormals.resize(E.size());

    for (int v = 0; v < dCN_V.size(); v++) {
        Verts.row(v) = dCN_V[v].new_pos.transpose();
        weights[v] = dCN_V[v].w;
    }

    for (int e = 0; e < E.size(); e++) {
        int he_pos = E[e].he;
        int he_neg = HE[he_pos].twin;
        Edges[e] = std::array<int, 2>({HE[he_neg].dest, HE[he_pos].dest});
        // Positive
        if (std::isfinite(heDefData[he_pos].newFrame.tangent.norm()) && std::isfinite(heDefData[he_pos].newFrame.l)) {
            posEdgeTangents[e] = Utils::eigenToGLM(heDefData[he_pos].newFrame.l * heDefData[he_pos].newFrame.tangent);
        } else {
            if (!std::isfinite(heDefData[he_pos].newFrame.tangent.norm())) {
                std::cout << "Bad Tangent Vector on HE " << he_pos << std::endl;
            }
            if (!std::isfinite(heDefData[he_pos].newFrame.l)) {
                std::cout << "Bad Length on HE " << he_pos << std::endl;
            }
            posEdgeTangents[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
        if (std::isfinite(heDefData[he_pos].newFrame.binormal.norm()) && std::isfinite(heDefData[he_pos].newFrame.w)) {
            posEdgeBinormals[e] = Utils::eigenToGLM(heDefData[he_pos].newFrame.w * heDefData[he_pos].newFrame.binormal);
        } else {
            if (!std::isfinite(heDefData[he_pos].newFrame.binormal.norm())) {
                std::cout << "Bad Binormal Vector on HE " << he_pos << std::endl;
            }
            if (!std::isfinite(heDefData[he_pos].newFrame.w)) {
                std::cout << "Bad Width on HE " << he_pos << std::endl;
            }
            posEdgeBinormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }

        if (std::isfinite(heDefData[he_pos].newFrame.normal.norm()) && std::isfinite(heDefData[he_pos].newFrame.h)) {
            posEdgeNormals[e] = Utils::eigenToGLM(heDefData[he_pos].newFrame.h * heDefData[he_pos].newFrame.normal);
        } else {
            if (!std::isfinite(heDefData[he_pos].newFrame.normal.norm())) {
                std::cout << "Bad Normal Vector on HE " << he_pos << std::endl;
            }
            if (!std::isfinite(heDefData[he_pos].newFrame.h)) {
                std::cout << "Bad Height on HE " << he_pos << std::endl;
            }
            posEdgeNormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }


        // Negative
        if (std::isfinite(heDefData[he_neg].newFrame.tangent.norm()) && std::isfinite(heDefData[he_neg].newFrame.l)) {
            negEdgeTangents[e] = Utils::eigenToGLM(heDefData[he_neg].newFrame.l * heDefData[he_neg].newFrame.tangent);
        } else {
            if (!std::isfinite(heDefData[he_neg].newFrame.binormal.norm())) {
                std::cout << "Bad Tangent Vector on HE " << he_neg << std::endl;
            }
            if (!std::isfinite(heDefData[he_neg].newFrame.w)) {
                std::cout << "Bad Length on HE " << he_neg << std::endl;
            }
            negEdgeTangents[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }

        if (std::isfinite(heDefData[he_neg].newFrame.binormal.norm()) && std::isfinite(heDefData[he_neg].newFrame.w)) {
            negEdgeBinormals[e] = Utils::eigenToGLM(heDefData[he_neg].newFrame.w * heDefData[he_neg].newFrame.binormal);
        } else {
            if (!std::isfinite(heDefData[he_neg].newFrame.binormal.norm())) {
                std::cout << "Bad Binormal Vector on HE " << he_neg << std::endl;
            }
            if (!std::isfinite(heDefData[he_neg].newFrame.w)) {
                std::cout << "Bad Width on HE " << he_neg << std::endl;
            }
            negEdgeBinormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }

        if (std::isfinite(heDefData[he_neg].newFrame.normal.norm()) && std::isfinite(heDefData[he_neg].newFrame.h)) {
            negEdgeNormals[e] = Utils::eigenToGLM(heDefData[he_neg].newFrame.h * heDefData[he_neg].newFrame.normal);
        } else {
            if (!std::isfinite(heDefData[he_neg].newFrame.normal.norm())) {
                std::cout << "Bad Normal Vector on HE " << he_neg << std::endl;
            }
            if (std::isfinite(heDefData[he_neg].newFrame.h)) {
                std::cout << "Bad Height on HE " << he_neg << std::endl;
            }
            negEdgeNormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
    }
    return 1;
}

int profilemover::cutmeshPolyscopeFormat(Eigen::MatrixXd& Verts, std::vector<std::vector<int>>& Faces,
                            std::vector<glm::vec3>& VertN, std::vector<glm::vec3>& FaceN,
                            std::vector<glm::vec3>& cornerIdx,
                            std::vector<glm::vec3>& projVecs) const {
    if (CM.active_v < 3) {
        return -1;
    }
    Verts.resize(CM.active_v, 3);
    Faces.resize(CM.active_f);
    VertN.resize(CM.active_v);
    FaceN.resize(CM.active_f);
    cornerIdx.resize(CM.active_v);
    projVecs.resize(CM.active_v);
    std::map<int, int> vToPSV;
    int curr_v = 0;
    for (int v = 0; v < CM.V.size(); v++) {
        if (CM.V[v].active) {
            if (curr_v >= CM.active_v) {
                std::cout << "Counted faces does not equal active faces" << std::endl;
                return -1;
            }
            Verts.row(curr_v) = CM.V[v].pos;
            vToPSV[v] = curr_v;
            VertN[curr_v] = Utils::eigenToGLM(CM.V[v].n);
            if (CM.cutData[v].label == 1 || CM.cutData[v].label == 2) {
                int c_idx = CM.cutData[v].corner_idx;
                if (c_idx < 0) {
                    std::cout << "Cut vertex has a non-existent corner index" << std::endl;
                }
                Eigen::Vector3d edgeMidpoint = 0.5 * (CM.V[CM.HE[c_idx].dest].pos + CM.V[CM.HE[CM.HE[c_idx].twin].dest].pos);
                cornerIdx[curr_v] = Utils::eigenToGLM(edgeMidpoint - CM.V[v].pos);
                projVecs[curr_v] = Utils::eigenToGLM(CM.cutData[v].defData.projVector);
            } else {
                cornerIdx[curr_v] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
                projVecs[curr_v] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
            }

            curr_v++;
        }
    }
    int curr_f = 0;
    for (int f = 0; f < CM.F.size(); f++) {
        if (CM.F[f].active) {
            if (curr_f >= CM.active_f) {
                std::cout << "Counted faces does not equal active faces" << std::endl;
                return -1;
            }
            std::vector<int> adjV = CM.faceAdjVertIdxs(f);
            for (int v = 0; v < adjV.size(); v++) {
                auto it = vToPSV.find(adjV[v]);
                if (it == vToPSV.end()) {
                    std::cout << "cutmeshPolyscopeFormat(): face references inactive/missing vertex "
                            << adjV[v] << " in face " << f << std::endl;
                    return -1;
                }
                adjV[v] = it->second;
            }
            Faces[curr_f] = adjV;
            FaceN[curr_f] = Utils::eigenToGLM(CM.F[f].n);
            curr_f++;
        }
    }
    return 1;
}

}   // namespace ProfileMover
