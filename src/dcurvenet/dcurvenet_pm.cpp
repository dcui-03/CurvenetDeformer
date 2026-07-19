#include "dcurvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <utility>
#include <iostream>

namespace DCurvenet {
/*
// Pre-compute maps
int dcurvenet::computedCNMaps(Eigen::SparseMatrix<double>& M_dCN_flat, 
                    Eigen::SparseMatrix<double>& M_3dCN_c,
                    const std::vector<std::pair<int, int>>& dCNV_to_c,
                    const std::vector<std::pair<int, int>>& dCNHE_to_c) {
    return 1;
}

// Compute dCN positions matrix
int dcurvenet::computedCNVerts(Eigen::MatrixXd& x_dCN, const std::vector<std::pair<int, int>>& dCNV_to_c) {
    #pragma omp parallel for
    for (int c = 0; c < dCNV_to_c.size(); c++) {
        x_dCN.row(dCNV_to_c[c].second) = V[dCNV_to_c[c].first].pos.transpose();
    }
    return 1;
}
// Compute the def grad operators at runtime
int dcurvenet::computeDefGradOperators(Eigen::MatrixXd& f_dCN_flat, Eigen::MatrixXd& f_dCN, 
                                       const std::vector<std::pair<int, int>>& dCNHE_to_c) {
    #pragma omp parallel for
    for (int c = 0; c < dCNHE_to_c.size(); c++) {
        const Eigen::Matrix3d& defGrad = HE[dCNHE_to_c[c].first].defData.defGrad;
        f_dCN_flat.row(dCNHE_to_c[c].second) = Utils::flattenMatrix3d(HE[dCNHE_to_c[c].first].defData.defGrad).transpose();
        f_dCN.block<3, 3>(3 * dCNHE_to_c[c].second, 0) = defGrad;
    }
    return 1;
}
*/

}   // namespace DCurvenet