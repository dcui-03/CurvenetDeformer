#include "utils.hpp"

#include <Eigen/Core>
#include <glm/vec3.hpp>
#include <vector>

namespace Utils {
    // GLM::vec3 to Eigen::Vector3d converter
    Eigen::Vector3d glmToEigen(const glm::vec3 input) {
        Eigen::Vector3d output;
        output(0) = static_cast<double>(input.x);
        output(1) = static_cast<double>(input.y);
        output(2) = static_cast<double>(input.z);
        return output;
    }

    // Eigen::Vector3d to GLM::vec3 converter
    glm::vec3 eigenToGLM(const Eigen::Vector3d input) {
        glm::vec3 output;
        output.x = static_cast<float>(input(0));
        output.y = static_cast<float>(input(1));
        output.z = static_cast<float>(input(2));
        return output;
    }

    // Entire mesh conversion routine Eigen to GLM
    void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM) {
        GLM.clear();
        GLM.resize(Eig.size());
        for (int v = 0; v < Eig.size(); v++) {
            GLM[v] = eigenToGLM(Eig[v]);
        }
        return;
    }

    // Entire mesh conversion routine GLM to Eigen
    void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM) {
        Eig.clear();
        Eig.resize(GLM.size());
        for (int v = 0; v < GLM.size(); v++) {
            Eig[v] = glmToEigen(GLM[v]);
        }
        return;
    }

    // Copy positions and connectivity into a copied container
    void copyPositions(const std::vector<Eigen::Vector3d>& V_old, std::vector<Eigen::Vector3d>& V_new) {
        V_new.clear();
        V_new.resize(V_old.size());
        for (int v = 0; v < V_old.size(); v++) {
            Eigen::Vector3d new_v = {V_old[v](0), V_old[v](1), V_old[v](2)};
            V_new[v] = new_v;
        }
        return;
    }

    void copyConnectivity(const std::vector<std::vector<int>>& T_old, std::vector<std::vector<int>>& T_new) {
        T_new.clear();
        T_new.resize(V_old.size());
        for (int f = 0; f < T_old.size(); f++) {
            std::vector<int> f_idxs;
            for (int v = 0; v < T_old[f].size(); v++) {
                f_idxs.push_back(T[f][v]);
            }
            V_new[v] = f_idxs;
        }
        return;
    }

} // namespace Utils