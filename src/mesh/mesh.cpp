#include "mesh.hpp"

#include "dcurvenet/dcurvenet.hpp"
#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <limits>
#include <utility>

// Mesh class functions for initialization

namespace Mesh {

// Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
mesh::mesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List) {
    if (!initHalfEdgeMesh(V_List, F_List)) {
        throw std::runtime_error("Failed to initialize halfedge mesh.");
    }
    computeFNormalsAreas();
    computeVNormalsAreas();
    computeHeightFuncs();
    computeMeanE();
    computeBBoxDiag();
    return;
}

// Add a discrete Curvenetwork Pointer
bool mesh::applyDiscreteCurvenet(DCurvenet::dcurvenet* dCurvenet) {
    dCN = dCurvenet;
    dCN_initialized = true;
    return true;
}

// Assign a discrete curvenet index to a halfedge
bool mesh::assignDCNtoHE(int he, const int dCN_idx, bool positive) {
    HE[he].dCN_idx = dCN_idx;
    HE[he].dCN_sign = positive;
    return true;
}

bool mesh::applyMeshRef(mesh* MRef) {
    M_ref = MRef;
    M_ref_initialized = true;
    return true;
}

// Initializes the half edge mesh (Verts, Edges, Faces, Halfedges) from a vertex and face list
bool mesh::initHalfEdgeMesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List) {
    clearMesh();
    // Guard against empty meshes
    if (V_List.empty()) {
        return false;
    }

    // 1. Initialize All Vertices
    V.resize(V_List.size());
    for (int v = 0; v < V_List.size(); v++) {
        V[v].pos = V_List[v];
    }
    active_v = V.size();

    // 2. Initialize faces, edges, and interior halfedges
    F.resize(F_List.size());
    // The number of face corners is a maximum on the number of edges needed
    int numFaceCorners = 0;
    for (const std::vector<int>& faceVerts : F_List) {
        if (faceVerts.size() < 3) { // Not a valid face
            return false;
        }
        numFaceCorners += faceVerts.size();
    }
    HE.reserve(2 * numFaceCorners);
    E.reserve(numFaceCorners);

    for (int f = 0; f < F_List.size(); f++) {
        const std::vector<int>& fVerts = F_List[f];
        const int fSize = static_cast<int>(fVerts.size());

        // Check if each vertex in the face has a valid index
        for (int i = 0; i < fSize; ++i) {
            int vi = fVerts[i];
            int vj = fVerts[(i + 1) % fSize];
            if (vi < 0 || vi >= V.size()) {
                return false;
            }
            if (vj < 0 || vj >= V.size()) {
                return false;
            }
            if (vi == vj) { // Degenerate face
                return false;
            }
        }
        F[f].verts = fVerts;

        // Temporary list of face HE's
        std::vector<int> faceHEs(fSize, -1);
        // Create interior halfedges for this face
        for (int i = 0; i < fSize; i++) {
            int vi = fVerts[i];
            int vj = fVerts[(i + 1) % fSize];
            // Keys for the vertex pair to HE map
            std::pair<int, int> dirKey = std::make_pair(vi, vj);
            std::pair<int, int> oppKey = std::make_pair(vj, vi);
            // Check that we don't already have this edge. If so, then there's a duplicate
            if (vertPairToHE.find(dirKey) != vertPairToHE.end()) {
                return false;
            }
            // Insert the new halfedge into the list and set its attributes
            int heIdx = HE.size();
            HE.emplace_back();

            HE[heIdx].dest = vj;
            HE[heIdx].face = f;

            faceHEs[i] = heIdx;
            vertPairToHE[dirKey] = heIdx;

            // Store one outgoing halfedge for vi
            if (V[vi].he == -1) {
                V[vi].he = heIdx;
            }

            // Check whether the opposite direction already exists
            auto oppIt = vertPairToHE.find(oppKey);
            if (oppIt != vertPairToHE.end()) {  // If so, then hook the two halfedges up
                int oppHE = oppIt->second;
                // If the opposite halfedge already has a twin, then more than 2 face touch the edge (non-manifold)
                if (HE[oppHE].twin != -1) {
                    return false;
                }
                // If we run into an invalid or non-existent edge, then something is also wrong
                int eIdx = HE[oppHE].edge;
                if (eIdx < 0 || eIdx >= E.size()) {
                    return false;
                }
                HE[heIdx].twin = oppHE;
                HE[oppHE].twin = heIdx;
                HE[heIdx].edge = eIdx;
            } else {    // If not, add in this new edge.
                int eIdx = E.size();
                E.emplace_back();

                E[eIdx].he = heIdx;
                HE[heIdx].edge = eIdx;
            }
        }

        // Wire next/prev HE's around the face.
        for (int i = 0; i < fSize; ++i) {
            int he = faceHEs[i];

            HE[he].next = faceHEs[(i + 1) % fSize];
            HE[he].prev = faceHEs[(i + fSize - 1) % fSize];
        }

        F[f].he = faceHEs[0];
    }

    // 3. Create boundary halfedges
    std::vector<int> boundaryHEs;
    const int numInteriorHEs = HE.size();
    // Iterate over interior halfedges and find any that have no twin (i.e., twin = -1)
    for (int he = 0; he < numInteriorHEs; ++he) {
        if (HE[he].twin != -1) {
            continue;
        }
        // Get the starting and ending vertices
        int u = HE[HE[he].prev].dest;
        int v = HE[he].dest;
        // Get an index and insert the new boundary half edge in + attributes
        int bhe = HE.size();
        HE.emplace_back();

        HE[bhe].dest = u;
        HE[bhe].twin = he;
        HE[bhe].edge = HE[he].edge;

        HE[he].twin = bhe;

        // Prefer boundary outgoing halfedge for boundary vertices (for easy querying)
        V[v].he = bhe;

        std::pair<int, int> bKey = std::make_pair(v, u);
        // If the boundary halfedge already exists somehow, then something is wrong
        if (vertPairToHE.find(bKey) != vertPairToHE.end()) {
            return false;
        }
        vertPairToHE[bKey] = bhe;
        boundaryHEs.push_back(bhe);
    }

    // 4. Connect next/prev for boundary halfedges
    std::map<int, int> boundaryOutgoingFromVertex;  // Store a map with the outgoing HE from each bdy vertex
    // Iterate over boundary halfedges
    for (int bhe : boundaryHEs) {
        // Boundary halfedge origin is the dest of its twin
        int origin = HE[HE[bhe].twin].dest;
        // If a vertex has more than one outgoing boundary halfedge, then it must be nonmanifold
        if (boundaryOutgoingFromVertex.find(origin) != boundaryOutgoingFromVertex.end()) {
            return false;
        }
        boundaryOutgoingFromVertex[origin] = bhe;
    }
    // Iterate over boundary vertices and find the associated next/prev
    for (int bhe : boundaryHEs) {
        int dest = HE[bhe].dest;
        // make sure that there exists an associated next vertex (i.e., valid boundary configuration)
        auto nextIt = boundaryOutgoingFromVertex.find(dest);
        if (nextIt == boundaryOutgoingFromVertex.end()) {
            return false;
        }
        int nextBHE = nextIt->second;

        HE[bhe].next = nextBHE;
        HE[nextBHE].prev = bhe;
    }

    active_e = E.size();
    active_f = F.size();

    return true;    // success!
}

// Clear all mesh attributes
bool mesh::clearMesh() {
    // Clear all lists
    V.clear();
    HE.clear();
    E.clear();
    F.clear();
    vertPairToHE.clear();
    // Reset number of vertices
    active_v = 0;
    active_e = 0;
    active_f = 0;
    // Reset mesh qualities
    meanE = 0.0;
    bboxDiag = 0.0;
    dCN = nullptr;
    return;
}

// Getters
Eigen::Vector3d& mesh::getVPos(int v) const {
    Eigen::Vector3d pos = V[v].pos;
    return pos;
}
Eigen::Vector3d& mesh::getVNormal(int v) const {
    Eigen::Vector3d n = V[v].n;
    return n;
}
Eigen::Vector3d& mesh::getFNormal(int f) const {
    Eigen::Vector3d n = F[f].n;
    return n;
}

// Get mean edge length
double mesh::getMeanE() const {
    return meanE;
}

 // Get bbox diagonal
double mesh::getBBoxDiag() const {
    return bboxDiag;
}

Eigen::VectorXd mesh::computeFaceHeight(int f) {
    const std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(f);
    int fSize = fVertsPos.size();
    Eigen::VectorXd faceH = Eigen::VectorXd::Zero(fSize);
    // Special handling for triangles (must be planar)
    if (fSize == 3) {
        faceH = Eigen::VectorXd({0.0, 0.0, 0.0});
        return faceH;
    }
    // 1. compute barycenter and face normal
    Eigen::Vector3d faceCenter = DECUtils::computeBarycenter(fVertsPos);
    Eigen::Vector3d faceN = F[f].n;
    // 2. Project face vertices onto the Newell plane and grab height
    vector2dList proj_v(fSize);
    // Build a basis
    Eigen::Vector3d t1;
    Eigen::Vector3d t2;
    Utils::buildPlaneBasis(faceN, t1, t2);
    for (int v = 0; v < fSize; v++) {
        Eigen::Vector3d proj3d = Utils::projectPointOntoPlane(faceN, faceCenter, fVertsPos[v]);
        // vector from old point to plane point
        Eigen::Vector3d heightVec = fVertsPos[v] - proj3d;
        // Get 2D version
        proj_v[v] = Utils::convertTo2D(proj3d, faceCenter, t1, t2);
        double height = heightVec.norm();   // How far we are from the plane
        if (height <= 1e-6) {   // If we are on/close to the surface, just snap to the plane
            faceH(v) = 0.0;
        } else if ((heightVec.normalized()).dot(faceN) > 0.0) { // We are above the plane
            faceH(v) = height;
        } else {    // We are below the plane
            faceH(v) = -1.0 * height;
        }
    }
    return faceH;
}

// Function which computes a single face's normal/area
double mesh::computeFVectorArea(int f, Eigen::Vector3d& fN) {
    std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(f);
    return DECUtils::vectorArea(fVertsPos, fN);
}

// Internal function to precompute normals on all mesh structures
void mesh::computeFNormalsAreas() {
    for (int f = 0; f < F.size(); f++) {
        if(!F[f].active) {
            continue;
        }
        Eigen::Vector3d fN = Eigen::Vector3d::Zero();
        F[f].fArea = computeFVectorArea(f, fN);
        F[f].n = fN;
    }
    return;
}

// Function which computes a single vertex's normal/area
double mesh::computeVNormalArea(int v, Eigen::Vector3d& vN, bool weight_fN) {
    // Iterate around the adjacent faces
    std::vector<int> fList = vertAdjFaces(v);
    vN = Eigen::Vector3d::Zero();
    double vArea = 0.0;
    // Iterate over face list and accumulate areas and normals
    for (int i = 0; i < fList.size(); i++) {
        int f = fList[f];
        double fArea = F[f].fArea/(F[f].verts.size());
        if (weight_fN) {
            vN += fArea * F[f].n;
        } else {
            vN += F[f].n;
        }
        vArea += fArea;
    }
    vN.normalize();
    return vArea;
}

void mesh::computeVNormalsAreas(bool weight_fN) {
    for (int v = 0; v < active_v; v++) {
        if (!V[v].active) {
            continue;
        }
        Eigen::Vector3d vN = Eigen::Vector3d::Zero();
        double vArea = computeVNormalArea(v, vN, weight_fN);
        V[v].n = vN;
        V[v].vArea = vArea;
    }
    return;
}

// Computes mean edge length on the mesh
void mesh::computeMeanE() {
    meanE = 0.0;
    if (active_e == 0) {
        return;
    }
    for (int e = 0; e < E.size(); ++e) {
        if (!E[e].active) {
            continue;
        }
        std::pair<int, int> eVerts = edgeAdjVerts(e);
        meanE += (V[eVerts.first].pos - V[eVerts.second].pos).norm();
    }
    meanE /= active_e;
    return;
}

// Computes the diagonal length of the mesh's AABB
void mesh::computeBBoxDiag() {
    int start = 0;
    Eigen::Vector3d minV;
    Eigen::Vector3d maxV;
    for (int v = 0; v < V.size() - 1; v++) {
        start = v;
        if (V[v].active) {
            minV = V[v].pos;
            maxV = V[v].pos;
            break;
        }
    }
    // Note: 0 or 1 point will collapse the bbox. Return an error for debug
    if ((active_v == 0) || (start >= V.size() - 1)) {
        // std::cout << "No BBox computable. Too few active vertices." << std::endl;
        return;
    }
    // Find most extreme points in mesh
    for (int v = start + 1; v < V.size(); v++) {
        minV = minV.cwiseMin(V[v].pos);
        maxV = maxV.cwiseMax(V[v].pos);
    }
    // get norm of the most extreme points
    bboxDiag = (maxV - minV).norm();
    return;
}

}   // namespace Mesh