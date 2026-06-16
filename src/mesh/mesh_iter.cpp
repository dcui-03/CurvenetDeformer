// File for mesh iterators and querying
#include "mesh.hpp"

#include <Eigen/Core>
#include <vector>
#include <utility>


namespace Mesh {

// Returns a CCW list of a vertex's OUTGOING halfedge indices
std::vector<int> mesh::vertAdjHEs(int v) const {
    std::vector<int> outgoingHEs;
    const int he0 = V[v].he;
    int he_curr = he0;
    do {
        outgoingHEs.push_back(he_curr);
        he_curr = HE[HE[he_curr].prev].twin;
    } while (he_curr != he0);
    return outgoingHEs;
}

// Returns a CCW list of a vertex's adjacent vertices
std::vector<int> mesh::vertAdjVerts(int v) const {
    std::vector<int> adjHE = vertAdjHEs(v);
    std::vector<int> adjVerts(adjHE.size());
    for (int he = 0; he < adjHE.size(); he++) {
        adjVerts[he] = HE[adjHE[he]].dest;
    }
    return adjVerts;
}

// Returns a CCW list of a vertex's adjacent faces
// INCLUDES BOUNDARY if a boundary is adjacent
std::vector<int> mesh::vertAdjFaces(int v) const {
    std::vector<int> adjHE = vertAdjHEs(v);
    std::vector<int> adjFaces(adjHE.size());
    for (int he = 0; he < adjHE.size(); he++) {
        int f = HE[adjHE[he]].face;
        bool faceExists = false;
        // Check that we have not already recorded this face
        for (int fi = 0; fi < adjFaces.size(); fi++) {
            if (f == adjFaces[fi]) {
                faceExists = true;
                break;
            }
        }
        if (faceExists) {   // If so, skip this case.
            continue;
        }
        adjFaces.push_back(f);
    }
    return adjFaces;
}

// Get an arbitrary vertex loop starting from a certain halfedge index
std::vector<int> mesh::vertLoop(int he) const {
    int max_search = active_e/2;
    int iter = 0;
    std::vector<int> vLoop;
    int he_curr = he;
    do {
        vLoop.push_back(HE[he_curr].dest);
        he_curr = HE[he_curr].next;
        iter++;
    } while((he_curr != he) && (iter <= max_search));
    return vLoop;
}

// Returns the endpoints of an edge in an arbitrary order.
std::pair<int, int> mesh::edgeAdjVerts(int e) const {
    int he = E[e].he;
    std::pair<int, int> v_pair = std::make_pair(HE[HE[he].prev].dest, HE[he].dest);
    return v_pair;
}

// Returns the adjacent face(s) of an edge (-1 indicates boundary)
std::pair<int, int> mesh::edgeAdjFaces(int e) const {
    int he0 = E[e].he;
    int he1 = HE[he0].twin;
    std::pair<int, int> v_pair = std::make_pair(HE[he0].face, HE[he1].face);
    return v_pair;
}

// Returns the halfedge index given the face index and edge index
// If both halfedges point to the same face, returns an arbitrary one
int mesh::halfedgeAtFaceEdge(int f, int e) const {
    int he0 = E[e].he;
    int he1 = HE[he0].twin;
    if (HE[he0].face == f) {
        return he0;
    }
    return he1;
}

// Get an arbitrary halfedge loop starting from a certain index
std::vector<int> mesh::halfedgeLoop(int he) const {
    int max_search = active_e/2;
    int iter = 0;
    std::vector<int> heLoop;
    int he_curr = he;
    do {
        heLoop.push_back(he_curr);
        he_curr = HE[he_curr].next;
        iter++;
    } while((he_curr != he) && (iter <= max_search));
    
    return heLoop;
}

// Returns a CCW list of a face's vertices
// NOTE: In case of scrambled vertex ordering (ex. interior loops), it's safest to do this by halfedge
std::vector<Eigen::Vector3d> mesh::faceAdjVerts(int f) const {
    std::vector<int> fVerts = faceAdjVertIdxs(f);
    std::vector<Eigen::Vector3d> fVertsPos(fVerts.size());
    for (int v = 0; v < fVerts.size(); v++) {
        fVertsPos[v] = V[fVerts[v]].pos;
    }
    return fVertsPos;
}
// Overload if given an fVerts
std::vector<Eigen::Vector3d> mesh::faceAdjVerts(std::vector<int> fVerts) const {
    std::vector<Eigen::Vector3d> fVertsPos(fVerts.size());
    for (int v = 0; v < fVerts.size(); v++) {
        fVertsPos[v] = V[fVerts[v]].pos;
    }
    return fVertsPos;
}

// Returns a CCW list of face vertex indices
std::vector<int> mesh::faceAdjVertIdxs(int f) const {
    std::vector<int> fHalfEdges = faceAdjHalfEdges(f);
    std::vector<int> fVerts(fVerts.size());
    for (int he = 0; he < fHalfEdges.size(); he++) {
        fVerts[he] = HE[fHalfEdges[he]].dest;
    }
    return fVerts;
}

// Returns a CCW list of a face's halfedges
std::vector<int> mesh::faceAdjHalfEdges(int f) const {
    std::vector<int> fHalfEdges;
    const int he0 = F[f].he;
    int he_curr = he0;
    do {
        fHalfEdges.push_back(he_curr);
        he_curr = HE[he_curr].next;
    } while(he_curr != he0);
    return fHalfEdges;
}

// Returns the outgoing boundary HE if a vertex is a boundary vertex, else returns -1
int mesh::vertIsBoundary(int v, bool fast) const {
    // Fast check: grab the outgoing halfedge
    if (fast) {
        if (HE[V[v].he].face == -1) {
            return true;
        } else {
            return false;
        }
    }
    // Hard check: Check all outgoing halfedges
    // Safety in case the soft boundary halfedge rule is accidentally violated
    const int he0 = V[v].he;
    int he_curr = he0;
    do {
        if (HE[he_curr].face == -1) {
            return true;
        }
        he_curr = HE[HE[he0].prev].twin;
    } while (he_curr != he0);
    return false;
}
// Returns true if halfedge is on boundary
bool mesh::halfedgeIsBoundary(int he) const {
    // Super simple face check
    if (HE[he].face == -1) {
        return true;
    }
    return false;
}

}   // namespace Mesh