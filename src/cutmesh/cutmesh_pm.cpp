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

int cutmesh::computeCMatrix(Eigen::SparseMatrix<double>& C_mat, std::vector<int>& cToCM, std::map<int, std::vector<int>>& mToC, const std::vector<int>& heToCMhe) {
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
    std::vector<T> tripletList;
    tripletList.reserve(num_he);
    
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
// TODO: Since each operates on a separate row of defGrads, is the parallelism safe?
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
            return -1;
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
        Eigen::Vector3d proj = V[v].defData.projVector;
        Eigen::Vector3d new_pos = -1 * V[v].defData.defGrad * V[v].defData.projVector;
        int dCN_corner = HE[V[v].corner_idx].dCN_idx;
        if (dCN_corner < 0) {
            std::cout << "Incorrect corner index assignment found in cutmesh::estimateCNPositions" << std::endl;
            return -1;
        }
        // Compute the vertex deformation gradient based on the corresponding dCN def grad
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
// TODO: Can we parallelize? If so, how?
// Problem is, we are trying to modify various rows of the deformed face final matrix
// They shouldn't collide since they're halfedges, but still... check if safe.
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