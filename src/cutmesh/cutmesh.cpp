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

// Cutmesh class functions for initialization

namespace Mesh {

// Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
cutmesh::cutmesh(mesh* MRef, DCurvenet::dcurvenet* dCN): M(MRef), dCN(dCN) {
    // Copy in vertex and edge data from the reference mesh
    if (!copyFromMesh()) {
        throw std::runtime_error("Failed to initialize copy from mesh.");
    }
    // Embed the dCN curves into this mesh
    if (embedCurves() != 1) {
        throw std::runtime_error("Failed to embed curves.");
    }
    // Actually apply cuts to the mesh
    /*
    if (cutMesh() != 1) {
        throw std::runtime_error("Failed to cut mesh.");
    }
    */
    return;
}

cutmesh::cutmesh() {

}

// Add a discrete Curvenetwork Pointer
bool cutmesh::applyDiscreteCurvenet(DCurvenet::dcurvenet* dCurvenet) {
    dCN = dCurvenet;
    return true;
}

// Add a discrete Curvenetwork Pointer
bool cutmesh::applyMeshRef(mesh* MRef) {
    M = MRef;
    return true;
}

// Assign a discrete curvenet index to a halfedge
bool cutmesh::assignDCNtoHE(int he, const int dCN_idx) {
    if (he < 0 || he >= HE.size()) {
        return false;
    }
    HE[he].dCN_idx = dCN_idx;
    return true;
}

// Count number of active elements
void cutmesh::countNumActive() {
    active_v = 0;
    active_e = 0;
    active_f = 0;
    for (int f = 0; f < F.size(); f++) {
        if (F[f].active) {
            active_f++;
        }
    }

    for (int e = 0; e < E.size(); e++) {
        if (E[e].active) {
            active_e++;
        }
    }
    
    for (int v = 0; v < V.size(); v++) {
        if (V[v].active) {
            active_v++;
        }
    }
    return;
}

// Initializes by copying the vertex, edge, and halfedge data from the reference mesh
bool cutmesh::copyFromMesh() {
    clearMesh();
    if (M->F.size() < 1 || M->V.size() < 3 || M->E.size() < 3) {
        return false;
    }
    // Because the source mesh must have fixed topology + no inactive attributes, we can copy directly without worrying about indexing issues
    V = M->V;
    HE = M->HE;
    E = M->E;

    // You can either copy faces or clear them.
    // Since you later rebuild faces, clearing is fine.
    F.clear();

    // If you want to intentionally invalidate old face IDs:
    for (int he = 0; he < HE.size(); he++) {
        HE[he].face = -1;
    }

    vertPairToHE = M->vertPairToHE;

    meanE = M->meanE;
    bboxDiag = M->bboxDiag;

    active_v = M->active_v;
    active_e = M->active_e;
    active_f = 0;

    // We will need to re-init faces later anyways, so we can just ignore for now
    return true;
}

}   // namespace Mesh