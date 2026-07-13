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

// Mesh class functions for various stages of the profile mover algorithm

namespace Mesh {


// Compute a mapping from halfedge indices in V and C to the cutmesh's halfedge indices
void cutmesh::computeHEMap(std::vector<int>& heToCMhe, std::map<int, int>& CMheTohe) {
    heToCMhe.clear();
    CMheTohe.clear();
    int num_he = 0;
    for (int f = 0; f < F.size(); f++) {
        std::vector<int> adjHE = faceAdjHalfEdges(f);
        for (int he = 0; he < adjHE.size(); he++) {
            if (HE[adjHE[he]].active && !HE[adjHE[he]].boundary) {
                heToCMhe.push_back(adjHE[he]);
                CMheTohe[adjHE[he]] = num_he;
                num_he++;
            }
        }
    }
    return;
}

// Compute a matrix V and also a mapping from V indices to the cutmesh cutvert indices
// Also fills in the mapping from mesh vert to V.
Eigen::SparseMatrix<double> cutmesh::computeVMatrix(std::vector<int>& vToCM, std::map<int, int>& mToV, const std::vector<int>& heToCMhe) {
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
    Eigen::SparseMatrix<double> V_mat(num_he, num_v);
    std::vector<T> tripletList;
    tripletList.reserve(num_he);
    
    // Fill the matrix
    // Our convention is that each halfedge uses its origin vertex
    for (int he = 0; he < heToCMhe.size(); he++) {
        int dest = HE[HE[heToCMhe[he]].twin].dest;
        if (V[dest].label == 0) {   // Get V only
            tripletList.push_back(T(heToCMhe[he], CMtoV[dest], 1.0));
        }
    }

    V_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return V_mat;
}

Eigen::SparseMatrix<double> cutmesh::computeCMatrix(std::vector<int>& cToCM, std::map<int, std::vector<int>>& mToC, const std::vector<int>& heToCMhe) {
    // Def. triplets for filling sparse matrices
    typedef Eigen::Triplet<double> T;
    cToCM.clear();
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
    Eigen::SparseMatrix<double> C_mat(num_he, num_c);
    std::vector<T> tripletList;
    tripletList.reserve(num_he);
    
    // Fill the matrix
    // Our convention is that each halfedge uses its origin vertex
    for (int he = 0; he < heToCMhe.size(); he++) {
        int dest = HE[HE[heToCMhe[he]].twin].dest;
        if (V[dest].label != 0) {
            tripletList.push_back(T(heToCMhe[he], CMtoC[dest], 1.0));
        }
    }

    C_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return C_mat;
}

// Compute the halfedge-based laplacian 
Eigen::SparseMatrix<double> cutmesh::computeHELaplacian(std::map<int, int>& CMheTohe) {
    typedef Eigen::Triplet<double> T;
    Eigen::SparseMatrix<double> L;
    std::vector<T> tripletList;
    tripletList.reserve(CMheTohe.size());
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
            int f_row = CMheTohe[adjHE_idxs[i]];
            for (int j = 0; j < faceL.cols(); j++) {
                int f_col = CMheTohe[adjHE_idxs[j]];
                tripletList.push_back(T(f_row, f_col, faceL(i, j)));
            }
        }
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
    return L;
}

// Compute deformation gradients on cut-vertices
// TODO: Since each operates on a separate row of defGrads, is the parallelism safe?
Eigen::MatrixXd cutmesh::computeDefGrads(const std::vector<int>& cToCM) {
    Eigen::MatrixXd defGrads(cToCM.size(), 9);
    // Grab deformation gradients from the corresponding cutmesh
    #pragma omp parallel for
    for (int c = 0; c < cToCM.size(); c++) {
        int v = cToCM[c];
        int dCN_prev = V[v].corner_idx;
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
    return defGrads;
}

// Apply solved deformation gradients to the cutmesh
void cutmesh::applyDefGrads(Eigen::MatrixXd defGrads, const std::vector<int>& vToCM) {
    #pragma omp parallel for
    for (int v = 0; v < vToCM.size(); v++) {
        Eigen::Matrix3d defGrad = Utils::compressVector9d(defGrads.row(v).transpose());
        V[vToCM[v]].defData.defGrad = defGrad;
    }
    return;
}

// Estimate projected curvenet positions
Eigen::MatrixXd cutmesh::estimateCNPositions(const std::vector<int>& cToCM) {
    Eigen::MatrixXd cnPos(cToCM.size(), 3);
    std::vector<Eigen::Vector3d> cnPos_vector(cToCM.size());

    #pragma omp parallel for
    // Grab deformation gradients from the corresponding cutmesh
    for (int c = 0; c < cToCM.size(); c++) {
        int v = cToCM[c];
        Eigen::Vector3d proj = V[v].defData.projVector;
        Eigen::Vector3d new_pos = -1 * V[v].defData.defGrad * V[v].defData.projVector;
        int dCN_corner = V[v].corner_idx;
        // Compute the vertex deformation gradient based on the corresponding dCN def grad
        if (V[v].label == 1) {
            Eigen::Vector3d target_pos = dCN->V[dCN->HE[dCN_corner].dest].pos;
            new_pos += target_pos;
        } else if (V[v].label == 2) {
            Eigen::Vector3d next_pos = dCN->V[dCN->HE[dCN_corner].dest].pos;
            Eigen::Vector3d prev_pos = dCN->V[dCN->HE[dCN->HE[dCN_corner].twin].dest].pos;
            new_pos += 0.5 * (next_pos + prev_pos);
        }
        cnPos_vector[c] = new_pos.transpose();
    }

    for (int c = 0; c < cToCM.size(); c++) {
        cnPos.row(c) = cnPos_vector[c];
    }
    return cnPos;
}

// Estimate the deformed faces
// TODO: Can we parallelize? If so, how?
// Problem is, we are trying to modify various rows of the deformed face final matrix
// They shouldn't collide since they're halfedges, but still... check if safe.
Eigen::MatrixXd cutmesh::estimateFaceDeformations(const std::map<int, int>& CMheTohe) {
    Eigen::MatrixXd deformedFaces(CMheTohe.size(), 3);

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
        Eigen::MatrixXd defFace = V_Matrix * defGrad.transpose();

        for (int v = 0; v < adjV.size(); v++) {
            int he = adjHE_idxs[v];
            deformedFaces.row(CMheTohe.at(he)) = defFace.row(v);
        }
    }
    return deformedFaces;
}

}   // namespace Mesh