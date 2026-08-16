#include "cutmesh.hpp"

// Cutmesh class functions for initialization

namespace Mesh {

// Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
cutmesh::cutmesh(mesh* MRef, Polynet::dcurvenet* dCN): M(MRef), dCN(dCN) {
    // Copy in vertex and edge data from the reference mesh
    if (!copyFromMesh()) {
        throw std::runtime_error("Failed to initialize copy from mesh.");
    }
    // Embed the dCN curves into this mesh
    if (embedCurves() != 1) {
        throw std::runtime_error("Failed to embed curves.");
    }
    // Actually apply cuts to the mesh
    if (cutMesh() != 1) {
        throw std::runtime_error("Failed to cut mesh.");
    }
    return;
}

cutmesh::cutmesh() {

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

    // Seed cutData for the copied-in original mesh vertices (each maps to itself)
    cutData.assign(V.size(), CutData());
    for (int v = 0; v < V.size(); v++) {
        cutData[v].projData.elType = 0;
        cutData[v].projData.elIdx = v;
    }

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