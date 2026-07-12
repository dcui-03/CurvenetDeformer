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

// Mesh class functions for initialization

namespace Mesh {

// Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
cutmesh::cutmesh(mesh* MRef, DCurvenet::dcurvenet* dCN): M(MRef), dCN(dCN) {
    // Copy in vertex and edge data from the reference mesh
    if (!copyFromMesh()) {
        throw std::runtime_error("Failed to initialize copy from mesh.");
    }
    // Embed the dCN curves into this mesh
    embedCurves();
    // Actually apply cuts to the mesh
    cutMesh();
    return;
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
    // Iterate over vertices and copy in their data
    for (int v = 0; v < M->V.size(); v++) {
        V.emplace_back();
        // TODO
    }

    // Iterate over halfedges and copy in their data
    for (int he = 0; he < M->HE.size(); he++) {
        // TODO
    }

    // Iterate over edges and copy in their data
    for (int e = 0; e < M->E.size(); e++) {
        // TODO
    }

    // We will need to re-init faces later anyways, so we can just ignore for now
    return true;
}

}   // namespace Mesh