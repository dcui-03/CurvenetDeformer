#include "cutmesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <limits>
#include <utility>
#include <algorithm>

// Utility functions for mesh (projection, insertion, etc.)

namespace Mesh {

// Inserts a vertex at a location and face
// Returns the index of the new vertex
// NOTE: The new vertex has no normal information or associated halfedge.
int cutmesh::insertVertex(Eigen::Vector3d pos, Eigen::Vector3d n, int label, int cornerIdx, int ref_Type, int ref_Idx, Eigen::Vector3d proj, Eigen::Matrix3d defGrad) {
    Vert newV = createVertex(pos, n, label, cornerIdx, ref_Type, ref_Idx, proj, defGrad);
    int v = V.size();
    V.push_back(newV);
    return v;
}

int cutmesh::insertVertex(Eigen::Vector3d pos, Eigen::Vector3d n, int label, int cornerIdx, vertProjData projData, vertDeformData defData) {
    Vert newV = createVertex(pos, n, label, cornerIdx, projData, defData);
    int v = V.size();
    V.push_back(newV);
    return v;
}

int cutmesh::insertVertex(Vert splitV) {
    int v = V.size();
    V.push_back(splitV);
    return v;
}

// Topologically splits an existing edge by adding a new vertex.
// NOTE: Added vertex does NOT need to lie on the edge
int cutmesh::splitEdge(int e, int new_v) {
    if (e >= E.size() || !E[e].active) {
        return -1;
    }
    if (new_v >= V.size() || !V[new_v].active) {
        return -1;
    }

    // Make and modify a copy of the current halfedges
    int he0_idx = E[e].he;
    int he1_idx = HE[he0_idx].twin;
    HalfEdge he0 = HE[he0_idx];
    HalfEdge he1 = HE[he1_idx];
    int u = HE[he1_idx].dest;
    int v = HE[he0_idx].dest;

    // Insert the new halfedges
    int he0_new = HE.size();
    int he1_new = he0_new + 1;
    HE.emplace_back();
    HE.emplace_back();
    HE[he0_new] = he0;
    HE[he1_new] = he1;
    // Insert the new edge
    int new_e = E.size();
    E.emplace_back();
    E[new_e].he = he0_new;
    active_e++;

    // Re-wire the halfedges
    // NOTE: Since it's copied from he0 and he1, we don't need to rewire halfedge faces
    HE[he0_new].twin = he1_new;
    HE[he1_new].twin = he0_new;
    HE[he0_new].edge = new_e;
    HE[he1_new].edge = new_e;
    HE[he0_idx].dest = new_v;
    HE[he1_new].dest = new_v;
    HE[he0_new].dest = v;
    HE[he0_new].boundary = HE[he0_idx].boundary;
    HE[he1_new].boundary = HE[he1_idx].boundary;
    V[new_v].he = he0_new;

    HE[HE[he0_idx].next].prev = he0_new;
    HE[HE[he1_idx].prev].next = he1_new;
    HE[he0_idx].next = he0_new;
    HE[he0_new].prev = he0_idx;
    HE[he1_new].next = he1_idx;
    HE[he1_idx].prev = he1_new;

    // Modify vertPairToHE
    // Erase old, and insert replacements
    vertPairToHE.erase({u, v});
    vertPairToHE.erase({v, u});
    vertPairToHE[{u, new_v}] = he0_idx;
    vertPairToHE[{new_v, u}] = he1_idx;
    // New keys
    vertPairToHE[{new_v, v}] = he0_new;
    vertPairToHE[{v, new_v}] = he1_new;

    return new_e;
}

// Insert an edge into the mesh between two existing vertices
// Note, at this point mesh faces and halfedge Next/Prev are broken and should not be used
int cutmesh::insertEdge(int v0, int v1, int dCN_idx0, int dCN_idx1) {
    // Validate inputs
    if (v0 < 0 || v0 >= V.size() || !V[v0].active) {    // valid v0
        return -1;
    }
    if (v1 < 0 || v1 >= V.size() || !V[v1].active) {    // valid v1
        return -1;
    }
    if (v0 == v1) { // cannot draw edge from a point to itself
        return -1;
    }

    // Check if edge already exists
    if (vertPairToHE.find({v0, v1}) != vertPairToHE.end() || vertPairToHE.find({v1, v0}) != vertPairToHE.end()) {
        return -1;
    }

    // Create new edge between the vertices and wire
    int e = E.size();
    E.emplace_back();
    int he0 = HE.size();
    int he1 = he0+1;
    HE.emplace_back();
    HE.emplace_back();
    // Fill in attributes
    E[e].he = he0;
    HE[he0].edge = e;
    HE[he1].edge = e;
    HE[he0].twin = he1;
    HE[he1].twin = he0;
    HE[he0].dest = v1;
    HE[he1].dest = v0;
    HE[he0].dCN_idx = dCN_idx0;
    HE[he1].dCN_idx = dCN_idx1;
    V[v0].he = he1;
    V[v1].he = he0;

    // New keys
    vertPairToHE[{v0, v1}] = he0;
    vertPairToHE[{v1, v0}] = he1;

    return e;
}

}   // namespace Mesh