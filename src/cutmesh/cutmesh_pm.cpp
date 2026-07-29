#include "cutmesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <limits>
#include <utility>
#include <iostream>

// Mesh class functions for various stages of the profile mover algorithm

namespace Mesh {


// Compute a mapping from halfedge indices in V and C to the cutmesh's halfedge indices
void cutmesh::computeHEMap(std::vector<int>& heToCMhe, std::map<int, int>& CMheTohe) {
    heToCMhe.clear();
    CMheTohe.clear();
    int num_he = 0;
    for (int f = 0; f < F.size(); f++) {
        if (!F[f].active) {
            continue;
        }
        std::vector<int> adjHE = faceAdjHalfEdges(f);
        for (int he = 0; he < adjHE.size(); he++) {
            if (HE[adjHE[he]].active && !HE[adjHE[he]].boundary) {
                heToCMhe.push_back(adjHE[he]);
                CMheTohe[adjHE[he]] = num_he;
                num_he++;
            }
        }
    }
    //std::cout << "Total used halfedges: " << num_he << " of " << HE.size() << std::endl;
    return;
}

// Compute a matrix V and also a mapping from V indices to the cutmesh cutvert indices
// Also fills in the mapping from mesh vert to V.
int cutmesh::computeVMatrix(Eigen::SparseMatrix<double>& V_mat, std::vector<int>& vToCM, std::map<int, int>& mToV, const std::vector<int>& heToCMhe) {
    // Def. triplets for filling sparse matrices
    typedef Eigen::Triplet<double> T;
    vToCM.clear();
    mToV.clear();
    std::map<int, int> CMtoV;
    int num_v = 0;
    // First just get number of vertices
    for (int v = 0; v < V.size(); v++) {
        if (V[v].label == 0 && V[v].active) {
            vToCM.push_back(v);
            CMtoV[v] = num_v;
            mToV[V[v].projData.elIdx] = num_v;
            num_v++;
        }
    }
    // Get number of halfedges
    int num_he = heToCMhe.size();
    V_mat.resize(num_he, num_v);
    V_mat.setZero();
    std::vector<T> tripletList;
    tripletList.reserve(num_he);
    
    // Fill the matrix
    // Our convention is that each halfedge uses its origin vertex
    for (int he = 0; he < heToCMhe.size(); he++) {
        int origin = HE[HE[heToCMhe[he]].twin].dest;
        if (V[origin].label == 0) {   // Get V only
            tripletList.push_back(T(he, CMtoV[origin], 1.0));
        }
    }

    // std::cout << "V_mat has " << tripletList.size() << " non-zero entries" << std::endl;
    V_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return 1;
}

int cutmesh::computeCMatrix(Eigen::SparseMatrix<double>& C_mat, Eigen::MatrixXd& c_rest, std::vector<int>& cToCM, std::map<int, std::vector<int>>& mToC, const std::vector<int>& heToCMhe) {
    // Def. triplets for filling sparse matrices
    typedef Eigen::Triplet<double> T;
    cToCM.clear();
    mToC.clear();
    std::map<int, int> CMtoC;
    int num_c = 0;
    // First just get number of vertices
    for (int v = 0; v < V.size(); v++) {
        if (V[v].label != 0 && V[v].active) {
            cToCM.push_back(v);
            CMtoC[v] = num_c;
            if (V[v].projData.elType == 0) {
                mToC[V[v].projData.elIdx].push_back(num_c);
            }
            num_c++;
        }
    }
    // Get number of halfedges
    int num_he = heToCMhe.size();
    C_mat.resize(num_he, num_c);
    C_mat.setZero();
    c_rest.resize(num_c, 3);
    std::vector<T> tripletList;
    tripletList.reserve(num_he);

    for (int c = 0; c < cToCM.size(); c++) {
        int he_idx = V[cToCM[c]].corner_idx;
        int dCN_v_idx = dCN->HE[HE[he_idx].dCN_idx].dest;
        c_rest.row(c) = (dCN->V[dCN_v_idx].pos).transpose();
    }
    
    // Fill the matrix
    // Our convention is that each halfedge uses its origin vertex
    for (int he = 0; he < heToCMhe.size(); he++) {
        int origin = HE[HE[heToCMhe[he]].twin].dest;
        if (V[origin].label != 0) {
            tripletList.push_back(T(he, CMtoC[origin], 1.0));
        }
    }

    // std::cout << "C_mat has " << tripletList.size() << " non-zero entries" << std::endl;
    C_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return 1;
}

int cutmesh::computeProjMatrix(Eigen::MatrixXd& projVecs, const std::vector<int>& cToCM) {
    projVecs.resize(cToCM.size(), 3);
    for (int c = 0; c < cToCM.size(); c++) {
        projVecs.row(c) = (V[cToCM[c]]).defData.projVector.transpose();
    }
    return 1;
}

int cutmesh::compute_dCNMaps(Eigen::SparseMatrix<double>& M_dCN_flat,
                             Eigen::SparseMatrix<double>& M_dCN_c,
                             const std::vector<int>& cToCM) {
    using T = Eigen::Triplet<double>;
    const int num_C   = static_cast<int>(cToCM.size());
    const int num_dCNhe = dCN->numHalfedges();
    const int num_dCNv  = dCN->numVerts();

    std::vector<T> grad_tripletList;
    std::vector<T> pos_tripletList;
    grad_tripletList.reserve(2 * num_C);
    pos_tripletList.reserve(2 * num_C);

    for (int c = 0; c < num_C; c++) {
        int v = cToCM[c];
        int CM_he = V[v].corner_idx;
        int dCN_prev = HE[CM_he].dCN_idx;

        if (V[v].label == 1) {  // If we are on a label 1 cut-vert, take the average of the adjacent halfedges
            int dCN_next = dCN->HE[dCN_prev].next;
            int CM_he = V[v].corner_idx;
            // 0.5 * (prev defGrad + next defGrad)
            grad_tripletList.emplace_back(c, dCN_prev, 0.5);
            grad_tripletList.emplace_back(c, dCN_next, 0.5);
            // Get position
            int target_v = dCN->HE[dCN_prev].dest;
            pos_tripletList.emplace_back(c, target_v, 1.0);
        } else if (V[v].label == 2) {   // Otherwise, just take the parent halfedge
            grad_tripletList.emplace_back(c, dCN_prev, 1.0);
            // Get current position
            int dCN_twin = dCN->HE[dCN_prev].twin;
            int next_v = dCN->HE[dCN_prev].dest;
            int prev_v = dCN->HE[dCN_twin].dest;
            pos_tripletList.emplace_back(c, next_v, 0.5);
            pos_tripletList.emplace_back(c, prev_v, 0.5);
        }
    }
    // Set maps
    M_dCN_flat.resize(num_C, num_dCNhe);
    M_dCN_flat.setZero();
    M_dCN_flat.setFromTriplets(grad_tripletList.begin(), grad_tripletList.end());

    M_dCN_c.resize(num_C, num_dCNv);
    M_dCN_c.setZero();
    M_dCN_c.setFromTriplets(pos_tripletList.begin(), pos_tripletList.end());

    return 1;
}

int cutmesh::computeFaceOps(Eigen::SparseMatrix<double>& M_v_F,
                            Eigen::SparseMatrix<double>& M_c_F,
                            std::vector<int>& M_he_F,
                            Eigen::MatrixXd& x_h,
                            const std::vector<int>& vToCM,
                            const std::vector<int>& cToCM,
                            const std::vector<int>& heToCMhe,
                            const std::map<int, int>& CMheTohe) {
    using T = Eigen::Triplet<double>;
    const int num_he = heToCMhe.size();
    const int num_V = vToCM.size();
    const int num_C = cToCM.size();

    // Build temporary inverse maps
    std::vector<int> CMtoV(V.size(), -1);
    std::vector<int> CMtoC(V.size(), -1);
    for (int i = 0; i < num_V; i++) {
        CMtoV[vToCM[i]] = i;
    }
    for (int i = 0; i < num_C; i++) {
        CMtoC[cToCM[i]] = i;
    }

    // Raw active face index -> compact active face row
    std::vector<int> faceToRow(F.size(), -1);
    int num_F = 0;

    for (int f = 0; f < F.size(); f++) {
        if (!F[f].active) {
            continue;
        }
        faceToRow[f] = num_F++;
    }

    M_he_F.assign(num_he, -1);
    x_h.resize(num_he, 3);
    x_h.setZero();

    std::vector<T> v_tripletList;
    std::vector<T> c_tripletList;
    // Triangle/quads/etc. reserve estimate
    v_tripletList.reserve(3 * num_F);
    c_tripletList.reserve(3 * num_F);

    for (int f = 0; f < static_cast<int>(F.size()); f++) {
        if (!F[f].active) {
            continue;
        }
        int fRow = faceToRow[f];
        std::vector<int> adjHE = faceAdjHalfEdges(f);
        std::vector<int> adjV  = faceAdjVertIdxs(f, true); // origin order
        double w = 1.0 / static_cast<double>(adjV.size());
        // Average over the adjacent vertices
        for (int v = 0; v < adjV.size(); v++) {
            int CM_v = adjV[v];
            if (V[CM_v].label == 0) {
                int vRow = CMtoV[CM_v];
                v_tripletList.emplace_back(fRow, vRow, w);
            } else {
                int cRow = CMtoC[CM_v];
                c_tripletList.emplace_back(fRow, cRow, w);
            }
        }

        // Per-halfedge face row and rest origin position
        for (int v = 0; v < adjHE.size(); v++) {
            int rawHE = adjHE[v];
            auto it = CMheTohe.find(rawHE);
            if (it == CMheTohe.end()) {
                // Boundary/inactive halfedges are not part of the halfedge unknowns.
                continue;
            }
            int hRow = it->second;
            int origin = HE[HE[rawHE].twin].dest;
            M_he_F[hRow] = fRow;
            x_h.row(hRow) = V[origin].pos.transpose();
        }
    }
    // Sanity
    for (int he = 0; he < num_he; he++) {
        if (M_he_F[he] < 0) {
            std::cout << "computeFaceOps(): missing face row for halfedge row." << std::endl;
            return -1;
        }
    }

    M_v_F.resize(num_F, num_V);
    M_v_F.setZero();
    M_v_F.setFromTriplets(v_tripletList.begin(), v_tripletList.end());

    M_c_F.resize(num_F, num_C);
    M_c_F.setZero();
    M_c_F.setFromTriplets(c_tripletList.begin(), c_tripletList.end());

    return 1;
}

int cutmesh::computeAssemblyOps(Eigen::SparseMatrix<double>& M_v_M,
                                Eigen::SparseMatrix<double>& M_c_M,
                                const std::map<int, int>& mToV,
                                const std::map<int, std::vector<int>>& mToC,
                                int num_M, int num_V, int num_C) {
    using T = Eigen::Triplet<double>;

    std::vector<T> v_tripletList;
    std::vector<T> c_tripletList;
    v_tripletList.reserve(mToV.size());

    int cReserve = 0;
    for (const std::pair<int, std::vector<int>>& kv : mToC) {
        cReserve += kv.second.size();
    }
    c_tripletList.reserve(cReserve);

    for (const std::pair<int, int>& kv : mToV) {
        int m = kv.first;
        int vRow = kv.second;
        if (m < 0 || m >= num_M) {
            std::cout << "computeAssemblyOps(): invalid mToV entry." << std::endl;
            return -1;
        }
        v_tripletList.emplace_back(m, vRow, 1.0);
    }

    for (const std::pair<int, std::vector<int>>& kv : mToC) {
        int m = kv.first;
        const std::vector<int>& cRows = kv.second;
        if (m < 0 || m >= num_M) {
            std::cout << "computeAssemblyOps(): invalid mToC entry." << std::endl;
            return -1;
        }
        const double w = 1.0 / static_cast<double>(cRows.size());

        for (int cRow : cRows) {
            c_tripletList.emplace_back(m, cRow, w);
        }
    }

    M_v_M.resize(num_M, num_V);
    M_v_M.setZero();
    M_v_M.setFromTriplets(v_tripletList.begin(), v_tripletList.end());

    M_c_M.resize(num_M, num_C);
    M_c_M.setZero();
    M_c_M.setFromTriplets(c_tripletList.begin(), c_tripletList.end());
    return 1;
}

// Create a vector of weights per constraint vertex v
int cutmesh::computeWeightOps(Eigen::VectorXd& weights, const std::vector<int>& cToCM) {
    weights.resize(cToCM.size());
    for (int c = 0; c < cToCM.size(); c++) {
        int v = cToCM[c];
        if (!V[v].active || V[v].label == 0) {
            continue;
        }
        int dCN_prev = HE[V[v].corner_idx].dCN_idx;
        if (V[v].label == 1) {
            int dCN_v = dCN->HE[dCN_prev].dest;
            weights[c] = dCN->V[dCN_v].w;
        } else if (V[v].label == 2) {
            int dCN_twin = dCN->HE[dCN_prev].twin;
            int next_v = dCN->HE[dCN_prev].dest;
            int prev_v = dCN->HE[dCN_twin].dest;
            weights[c] = 0.5 * (dCN->V[next_v].w + dCN->V[prev_v].w);
        }
    }
    return 1;
}

// Compute the halfedge-based laplacian 
int cutmesh::computeHELaplacian(Eigen::SparseMatrix<double>& L, std::map<int, int>& CMheTohe) {
    typedef Eigen::Triplet<double> T;
    int num_he = CMheTohe.size();
    L.resize(num_he, num_he);
    L.setZero();
    int num_entries = 0;
    // Get an estimate of how many entries we will need
    for (int f = 0; f < F.size(); f++) {
        if (!F[f].active) {
            continue;
        }
        // Get the adjacent vertices and halfedges
        int n = faceAdjHalfEdges(f).size();
        num_entries += n*n;
    }
    std::vector<T> tripletList;
    tripletList.reserve(num_entries);
    // Compute face Laplacian
    for (int f = 0; f < F.size(); f++) {
        if (!F[f].active) {
            continue;
        }
        // Get the adjacent vertices and halfedges
        std::vector<int> adjHE_idxs = faceAdjHalfEdges(f);
        // Make sure to set the flag to origin so we can align easily with V and C's indexing
        std::vector<Eigen::Vector3d> adjV = faceAdjVerts(f, true);
        // Compute the face Laplacian
        Eigen::MatrixXd faceL = DECUtils::faceLaplacianOp(adjV);

        // Redistribute the Laplacian to its associated indices
        for (int i = 0; i < faceL.rows(); i++) {
            if (CMheTohe.find(adjHE_idxs[i]) == CMheTohe.end()) {
                std::cout << "Missing halfedge in CMheTohe during halfedge Laplacian assembly." << std::endl;
                return -1;
            }
            int f_row = CMheTohe.at(adjHE_idxs[i]);
            for (int j = 0; j < faceL.cols(); j++) {
                if (CMheTohe.find(adjHE_idxs[j]) == CMheTohe.end()) {
                    std::cout << "Missing halfedge in CMheTohe during halfedge Laplacian assembly." << std::endl;
                    return -1;
                }
                int f_col = CMheTohe.at(adjHE_idxs[j]);
                tripletList.push_back(T(f_row, f_col, faceL(i, j)));
            }
        }
    }
    // std::cout << "L has " << num_entries << " non-zero entries" << std::endl;
    L.setFromTriplets(tripletList.begin(), tripletList.end());
    return 1;
}

// Compute deformation gradients on cut-vertices
int cutmesh::computeDefGrads(Eigen::MatrixXd& defGrads, const std::vector<int>& cToCM) {
    // TODO: Assert so we don't have to resize
    defGrads.resize(cToCM.size(), 9);
    defGrads.setZero();
    // Grab deformation gradients from the corresponding cutmesh
    #pragma omp parallel for
    for (int c = 0; c < cToCM.size(); c++) {
        int v = cToCM[c];
        int dCN_prev = HE[V[v].corner_idx].dCN_idx;
        if (dCN_prev < 0) {
            std::cout << "Incorrect corner index assignment found in cutmesh::computeDefGrads" << std::endl;
            // return -1;
        }
        Eigen::Matrix3d defGrad;
        // Compute the vertex deformation gradient based on the corresponding dCN def grad
        if (V[v].label == 1) {
            int dCN_next = dCN->HE[dCN_prev].next;
            defGrad = 0.5 * (dCN->HE[dCN_prev].defData.defGrad + dCN->HE[dCN_next].defData.defGrad);
        } else if (V[v].label == 2) {
            defGrad = dCN->HE[dCN_prev].defData.defGrad;
        }
        V[v].defData.defGrad = defGrad;
        defGrads.row(c) = Utils::flattenMatrix3d(defGrad).transpose();
    }
    // std::cout << "Constraint def grads has size " << cToCM.size() << std::endl;
    return 1;
}

// Apply solved deformation gradients to the cutmesh
void cutmesh::applyDefGrads(const Eigen::MatrixXd& defGrads, const std::vector<int>& vToCM) {
    // std::cout << "Num def grads to distribute: " << defGrads.rows() << std::endl;
    #pragma omp parallel for
    for (int v = 0; v < vToCM.size(); v++) {
        Eigen::Matrix3d defGrad = Utils::compressVector9d(defGrads.row(v).transpose());
        V[vToCM[v]].defData.defGrad = defGrad;
    }
    // std::cout << "Distributed " << vToCM.size() << " def grads." << std::endl;
    return;
}

// Estimate projected curvenet positions
int cutmesh::estimateCNPositions(Eigen::MatrixXd& cnPos, const std::vector<int>& cToCM) {
    cnPos.resize(cToCM.size(), 3);
    cnPos.setZero();
    std::vector<Eigen::Vector3d> cnPos_vector(cToCM.size());

    #pragma omp parallel for
    // Grab deformation gradients from the corresponding cutmesh
    for (int c = 0; c < cToCM.size(); c++) {
        int v = cToCM[c];
        int dCN_corner = HE[V[v].corner_idx].dCN_idx;
        if (dCN_corner < 0) {
            std::cout << "Incorrect corner index assignment found in cutmesh::estimateCNPositions" << std::endl;
            // return -1;
        }
        Eigen::Vector3d new_pos = -1 * V[v].defData.defGrad * V[v].defData.projVector;
        // Compute the vertex deformation position based on the corresponding dCN def grad
        if (V[v].label == 1) {
            Eigen::Vector3d target_pos = dCN->V[dCN->HE[dCN_corner].dest].new_pos;
            new_pos += target_pos;
        } else if (V[v].label == 2) {
            Eigen::Vector3d next_pos = dCN->V[dCN->HE[dCN_corner].dest].new_pos;
            Eigen::Vector3d prev_pos = dCN->V[dCN->HE[dCN->HE[dCN_corner].twin].dest].new_pos;
            new_pos += 0.5 * (next_pos + prev_pos);
        }
        cnPos_vector[c] = new_pos.transpose();
    }

    for (int c = 0; c < cToCM.size(); c++) {
        cnPos.row(c) = cnPos_vector[c];
    }
    return 1;
}

// Estimate the deformed faces
int cutmesh::estimateFaceDeformations(Eigen::MatrixXd& deformedFaces, const std::map<int, int>& CMheTohe, bool arap) {
    deformedFaces.resize(CMheTohe.size(), 3);
    deformedFaces.setZero();

    // For deformed faces
    #pragma omp parallel for 
    for (int f = 0; f < F.size(); f++) {
        if (!F[f].active) {
            continue;
        }
        // Get the adjacent vertices and halfedges
        std::vector<int> adjHE_idxs = faceAdjHalfEdges(f);
        // Make sure to set the flag to origin so we can align easily with V and C's indexing
        std::vector<int> adjV_idxs = faceAdjVertIdxs(f, true);
        std::vector<Eigen::Vector3d> adjV = adjVerts(adjV_idxs);
        Eigen::MatrixXd V_Matrix = DECUtils::posOp(adjV);
        // Average the deformation gradient
        Eigen::Matrix3d defGrad = Eigen::Matrix3d::Zero();
        for (int v = 0; v < adjV.size(); v++) {
            defGrad += V[adjV_idxs[v]].defData.defGrad;
        }
        defGrad /= adjV.size();
        Eigen::MatrixXd defFace;
        // ARAP: Extract rotations (test)
        if (arap) {
            Eigen::Matrix3d R, S;
            Utils::polarDecomposition(defGrad, R, S);
            defFace = V_Matrix * R.transpose();
        } else {
            defFace = V_Matrix * defGrad.transpose();
        }

        for (int v = 0; v < adjV.size(); v++) {
            int he = adjHE_idxs[v];
            deformedFaces.row(CMheTohe.at(he)) = defFace.row(v);
        }
    }
    return 1;
}

}   // namespace Mesh