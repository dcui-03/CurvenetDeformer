#include "utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/point_mesh_squared_distance.h>
#include <glm/vec3.hpp>
#include <limits>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

namespace Utils {

// Helper for computing 
// NOTE: This is TEMPORARY, only for triangle meshes
int closestPointNormalOnMesh(
    const Eigen::Vector3d& p,
    const Eigen::MatrixXd& V,
    const std::vector<std::vector<int>>& faces,
    Eigen::Vector3d& n)
{
    n = Eigen::Vector3d::Zero();

    // Failure cases
    if (V.rows() == 0) {
        return -1;
    }
    for (const auto& f : faces) {
        if (f.size() != 3) {
            return -1;
        }
    }
    if (faces.empty()) {
        return -1;
    }

    // Convert std::vector<std::vector<int>> faces to Eigen::MatrixXi
    Eigen::MatrixXi F(faces.size(), 3);

    for (int i = 0; i < static_cast<int>(faces.size()); ++i) {
        for (int j = 0; j < 3; ++j) {
            int vid = faces[i][j];
            // Not explicitly requested, but prevents invalid memory access.
            if (vid < 0 || vid >= V.rows()) {
                return -1;
            }
            F(i, j) = vid;
        }
    }

    // Bounding box diagonal tolerance
    Eigen::Vector3d bbMin = V.colwise().minCoeff();
    Eigen::Vector3d bbMax = V.colwise().maxCoeff();

    double bboxDiag = (bbMax - bbMin).norm();
    double snapTol = 1e-6 * bboxDiag;

    // Query closest point on mesh
    Eigen::MatrixXd P(1, 3);
    P.row(0) = p.transpose();

    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;

    igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

    int closestFace = I[0];
    Eigen::Vector3d q = C.row(0).transpose();

    // Precompute face normals and double areas
    std::vector<Eigen::Vector3d> faceNormals(F.rows(), Eigen::Vector3d::Zero());
    std::vector<double> faceDoubleAreas(F.rows(), 0.0);

    for (int fi = 0; fi < F.rows(); ++fi) {
        Eigen::Vector3d a = V.row(F(fi, 0)).transpose();
        Eigen::Vector3d b = V.row(F(fi, 1)).transpose();
        Eigen::Vector3d c = V.row(F(fi, 2)).transpose();

        Eigen::Vector3d rawNormal = (b - a).cross(c - a);
        double doubleArea = rawNormal.norm();

        faceDoubleAreas[fi] = doubleArea;

        if (doubleArea > 0.0) {
            faceNormals[fi] = rawNormal / doubleArea;
        }
    }

    // Precompute area-weighted vertex normals
    std::vector<Eigen::Vector3d> vertexNormals(V.rows(), Eigen::Vector3d::Zero());

    for (int fi = 0; fi < F.rows(); ++fi) {
        Eigen::Vector3d areaWeightedNormal =
            faceDoubleAreas[fi] * faceNormals[fi];

        for (int lv = 0; lv < 3; ++lv) {
            vertexNormals[F(fi, lv)] += areaWeightedNormal;
        }
    }

    for (int vi = 0; vi < V.rows(); ++vi) {
        if (vertexNormals[vi].norm() > 0.0) {
            vertexNormals[vi].normalize();
        }
    }

    auto closestPointOnSegment = [](
        const Eigen::Vector3d& x,
        const Eigen::Vector3d& a,
        const Eigen::Vector3d& b) -> Eigen::Vector3d {
        Eigen::Vector3d ab = b - a;
        double denom = ab.squaredNorm();

        if (denom == 0.0) {
            return a;
        }

        double t = (x - a).dot(ab) / denom;
        t = std::max(0.0, std::min(1.0, t));

        return a + t * ab;
    };

    // First: snap to vertex if close enough
    for (int lv = 0; lv < 3; ++lv) {
        int vi = F(closestFace, lv);
        Eigen::Vector3d v = V.row(vi).transpose();

        if ((q - v).norm() <= snapTol) {
            n = vertexNormals[vi];

            if (n.norm() == 0.0) {
                n = faceNormals[closestFace];
            }

            return 1;
        }
    }

    // Second: snap to edge if close enough
    for (int le = 0; le < 3; ++le) {
        int v0 = F(closestFace, le);
        int v1 = F(closestFace, (le + 1) % 3);

        Eigen::Vector3d a = V.row(v0).transpose();
        Eigen::Vector3d b = V.row(v1).transpose();

        Eigen::Vector3d qEdge = closestPointOnSegment(q, a, b);

        if ((q - qEdge).norm() <= snapTol) {
            Eigen::Vector3d edgeNormal = Eigen::Vector3d::Zero();

            // Average normals of all faces adjacent to this edge.
            // For manifold meshes this is usually 1 or 2 faces.
            for (int fj = 0; fj < F.rows(); ++fj) {
                bool hasV0 = false;
                bool hasV1 = false;

                for (int k = 0; k < 3; ++k) {
                    if (F(fj, k) == v0) hasV0 = true;
                    if (F(fj, k) == v1) hasV1 = true;
                }

                if (hasV0 && hasV1) {
                    edgeNormal += faceNormals[fj];
                }
            }

            if (edgeNormal.norm() > 0.0) {
                n = edgeNormal.normalized();
            }
            else {
                n = faceNormals[closestFace];
            }
            return 1;
        }
    }

    // Otherwise, closest point is treated as being on the face interior
    n = faceNormals[closestFace];
    if (n.norm() == 0.0) {
        return -1;
    }
    return 1;
}

// HELPERS FOR CONVERSION/COPYING

// GLM::vec3 to Eigen::Vector3d converter
Eigen::Vector3d glmToEigen(const glm::vec3 input) {
    Eigen::Vector3d output;
    output(0) = static_cast<double>(input.x);
    output(1) = static_cast<double>(input.y);
    output(2) = static_cast<double>(input.z);
    return output;
}

// Eigen::Vector3d to GLM::vec3 converter
glm::vec3 eigenToGLM(const Eigen::Vector3d input) {
    glm::vec3 output;
    output.x = static_cast<float>(input(0));
    output.y = static_cast<float>(input(1));
    output.z = static_cast<float>(input(2));
    return output;
}

// Convert an eigen matrix with 3 columns to a std::vector
void EigM3toStdV(const Eigen::MatrixXd& mat, std::vector<Eigen::Vector3d>& vec) {
    // Assert that matrix columns = 3
    if (mat.cols() != 3) {
        return;
    }
    vec.resize(mat.rows());
    for (int v = 0; v < mat.rows(); v++) {
        vec[v] = mat.row(v).transpose();
    }
    return;
}

// Entire mesh conversion routine Eigen to GLM
void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM) {
    GLM.clear();
    GLM.resize(Eig.size());
    for (int v = 0; v < Eig.size(); v++) {
        GLM[v] = eigenToGLM(Eig[v]);
    }
    return;
}

// Entire mesh conversion routine GLM to Eigen
void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM) {
    Eig.clear();
    Eig.resize(GLM.size());
    for (int v = 0; v < GLM.size(); v++) {
        Eig[v] = glmToEigen(GLM[v]);
    }
    return;
}

// Copy positions and connectivity into a copied container
void copyPositions(const std::vector<Eigen::Vector3d>& V_old, std::vector<Eigen::Vector3d>& V_new) {
    V_new.clear();
    V_new.resize(V_old.size());
    for (int v = 0; v < V_old.size(); v++) {
        Eigen::Vector3d new_v = {V_old[v](0), V_old[v](1), V_old[v](2)};
        V_new[v] = new_v;
    }
    return;
}

void copyConnectivity(const std::vector<std::vector<int>>& T_old, std::vector<std::vector<int>>& T_new) {
    T_new.clear();
    T_new.resize(T_old.size());
    for (int f = 0; f < T_old.size(); f++) {
        std::vector<int> f_idxs;
        for (int v = 0; v < T_old[f].size(); v++) {
            f_idxs.push_back(T_old[f][v]);
        }
        T_new[f] = f_idxs;
    }
    return;
}

// SORTING
void doubleListIdxSort(std::vector<double>& ref_List, std::vector<int>& idx_List) {
    if (ref_List.size() != idx_List.size()) {
        return;
    }

    const int n = static_cast<int>(ref_List.size());

    for (int i = 0; i < n - 1; ++i) {
        bool swapped = false;

        for (int j = 0; j < n - i - 1; ++j) {
            if (ref_List[j] > ref_List[j + 1]) {
                std::swap(ref_List[j], ref_List[j + 1]);
                std::swap(idx_List[j], idx_List[j + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }
    }
}

// Insert an integer entry in a list between two specified values
bool insertIdxBetweenPair(std::vector<int>& idxList, int a, int b, int new_idx) {
    for (int i = 0; i < idxList.size(); ++i) {
        int j = (i + 1) % idxList.size();
        if (idxList[i] == a && idxList[j] == b) {
            idxList.insert(idxList.begin() + j, new_idx);
            return true;
        }
    }
    return false;
}

// Flattens an Eigen::Matrix3d into a 9x1 row vector
// NOTE: Does so column-wise!
Eigen::VectorXd flattenMatrix3d(const Eigen::Matrix3d& F) {
    return F.reshaped();
}

// Compresses a 9x1 Eigen::VectorXd into an Eigen::Matrix3d
// Assumes column-wise storage
Eigen::Matrix3d compressVector9d(const Eigen::VectorXd& f) {
    assert(f.size() == 9);
    return Eigen::Map<const Eigen::Matrix3d>(f.data());
}

// MESH HELPERS
std::pair<int, int> undirectedKey(int a, int b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

int edgeDirRelativeToKey(int a, int b) {
    return (a < b) ? +1 : -1;
}

bool orientFacesConsistently(std::vector<std::vector<int>>& F_List) {
    struct FaceEdgeUse {
        int face = -1;
        int localEdge = -1;
        int dir = 0; // +1 if stored as min->max, -1 if max->min
    };
    using EdgeKey = std::pair<int, int>;

    std::map<EdgeKey, std::vector<FaceEdgeUse>> edgeUses;

    // 1. Build undirected edge -> incident face uses.
    for (int f = 0; f < static_cast<int>(F_List.size()); f++) {
        const auto& face = F_List[f];
        int n = static_cast<int>(face.size());

        if (n < 3) {
            return false;
        }

        for (int i = 0; i < n; i++) {
            int a = face[i];
            int b = face[(i + 1) % n];

            if (a == b) {
                return false;
            }

            EdgeKey key = undirectedKey(a, b);
            edgeUses[key].push_back(FaceEdgeUse{
                f,
                i,
                edgeDirRelativeToKey(a, b)
            });
        }
    }

    // 2. Reject nonmanifold edges for this mesh structure.
    for (const auto& kv : edgeUses) {
        if (kv.second.size() > 2) {
            return false;
        }
    }

    // 3. Build face adjacency with "same direction?" relation.
    std::vector<std::vector<std::pair<int, bool>>> faceAdj(F_List.size());

    for (const auto& kv : edgeUses) {
        const auto& uses = kv.second;

        if (uses.size() != 2) {
            continue; // boundary edge
        }

        const FaceEdgeUse& a = uses[0];
        const FaceEdgeUse& b = uses[1];

        // If two faces use the shared undirected edge in the same direction,
        // one of them must be flipped.
        bool sameDir = (a.dir == b.dir);

        faceAdj[a.face].push_back({b.face, sameDir});
        faceAdj[b.face].push_back({a.face, sameDir});
    }

    // 4. BFS assign flip parity per connected component.
    // flip[f] == 0 means keep original orientation.
    // flip[f] == 1 means reverse this face.
    std::vector<int> flip(F_List.size(), -1);

    for (int root = 0; root < static_cast<int>(F_List.size()); root++) {
        if (flip[root] != -1) {
            continue;
        }

        flip[root] = 0;
        std::queue<int> q;
        q.push(root);

        while (!q.empty()) {
            int f = q.front();
            q.pop();

            for (auto [g, sameDir] : faceAdj[f]) {
                // If sameDir is true, neighbor must have opposite flip parity.
                // If sameDir is false, neighbor must have same flip parity.
                int requiredFlip = flip[f] ^ static_cast<int>(sameDir);

                if (flip[g] == -1) {
                    flip[g] = requiredFlip;
                    q.push(g);
                } else if (flip[g] != requiredFlip) {
                    // Contradiction: non-orientable or inconsistent connectivity.
                    return false;
                }
            }
        }
    }

    // 5. Apply flips.
    for (int f = 0; f < static_cast<int>(F_List.size()); f++) {
        if (flip[f]) {
            std::reverse(F_List[f].begin(), F_List[f].end());
        }
    }

    return true;
}

// GEOMETRY HELPERS

// Rotation-variant SVD
// Compute the Rotation Variant SVD
// Input an empty U, Sigma, V
void rotationVariantSVD(Eigen::Matrix3d& mat, Eigen::Matrix3d& U, Eigen::Vector3d& Sigma, Eigen::Matrix3d& V) {
    // Compute SVD to get Sigma, U and V
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    Sigma = svd.singularValues();
    V = svd.matrixV();

    // Compute L and remove reflections from U and V
    Eigen::Matrix3d L;
    L.setIdentity();
    L(2,2) = (U*V.transpose()).determinant();
    Sigma(2) = Sigma(2) * L(2,2); // To keep it a vector (taken from HOBAK)

    double u_det = U.determinant();
    double v_det = V.determinant();
    if (u_det < 0 && v_det > 0) {
        U = U * L;
    } else if (u_det > 0 && v_det < 0) {
        V = V * L;
    }

    return;
}

// Compute Polar Decomposition given rotation variant SVD
// Returns a vector containing R and then S
void polarDecomposition(Eigen::Matrix3d& mat, Eigen::Matrix3d& R, Eigen::Matrix3d& S) {
    Eigen::Matrix3d U, V; 
    Eigen::Vector3d Sigma;
    rotationVariantSVD(mat, U, Sigma, V);

    // Put together R and S from the inputs
    R = U * V.transpose();
    S = V * Sigma.asDiagonal() * V.transpose();
    return;
}

// Compute 3D signed angle between two vectors
double signedAngle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& axis, bool positive) {
    const double eps = 1e-12;

    if (v0.squaredNorm() <= eps || v1.squaredNorm() <= eps || axis.squaredNorm() <= eps) {
        return 0.0;
    }

    Eigen::Vector3d a = v0.normalized();
    Eigen::Vector3d b = v1.normalized();
    Eigen::Vector3d n = axis.normalized();

    double sinTheta = n.dot(a.cross(b));
    double cosTheta = std::clamp(a.dot(b), -1.0, 1.0);

    double sAngle = std::atan2(sinTheta, cosTheta);

    if (positive && sAngle < 0.0) {
        sAngle += 2.0 * M_PI;
    }

    return sAngle;
}

// Find the closest point to a triangle
Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d> triVerts, const Eigen::Vector3d p) {
    const double eps = 1e-8;
    const Eigen::Vector3d& a = triVerts[0];
    const Eigen::Vector3d& b = triVerts[1];
    const Eigen::Vector3d& c = triVerts[2];

    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);

    // Vertex region outside A
    if (d1 <= 0.0 && d2 <= 0.0) {
        return a;
    }

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);

    // Vertex region outside B
    if (d3 >= 0.0 && d4 <= d3) {
        return b;
    }

    // Edge region AB
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double denom = d1 - d3;
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, a, b, true);
        }

        const double t = d1 / denom;
        return a + t * ab;
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);

    // Vertex region outside C
    if (d6 >= 0.0 && d5 <= d6) {
        return c;
    }

    // Edge region AC
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double denom = d2 - d6;
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, a, c, true);
        }
        const double t = d2 / denom;
        return a + t * ac;
    }

    // Edge region BC
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double denom = (d4 - d3) + (d5 - d6);
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, b, c, true);
        }
        const double t = (d4 - d3) / denom;
        return b + t * (c - b);
    }

    // Inside face region
    const double denom = va + vb + vc;
    if (std::abs(denom) <= eps) {
        // Degenerate triangle fallback: closest point among three edges.
        Eigen::Vector3d pab = Utils::closestPointOnSegment3D(p, a, b, true);
        Eigen::Vector3d pac = Utils::closestPointOnSegment3D(p, a, c, true);
        Eigen::Vector3d pbc = Utils::closestPointOnSegment3D(p, b, c, true);

        double dab = (p - pab).squaredNorm();
        double dac = (p - pac).squaredNorm();
        double dbc = (p - pbc).squaredNorm();

        if (dab <= dac && dab <= dbc) return pab;
        if (dac <= dab && dac <= dbc) return pac;
        return pbc;
    }

    const double invDenom = 1.0 / denom;
    const double vBary = vb * invDenom;
    const double wBary = vc * invDenom;
    return a + vBary * ab + wBary * ac;
}

// Converts a 3D direction into an angle in the tangent plane spanned by t1, t2.
// Returns false if the projected direction is degenerate.
bool directionAngleInPlane(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& target,
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double& theta
) {
    const double eps = 1e-12;

    Eigen::Vector3d n = normal;
    if (n.squaredNorm() <= eps) {
        return false;
    }
    n.normalize();

    Eigen::Vector3d d = target - origin;
    d = d - d.dot(n) * n;

    if (d.squaredNorm() <= eps) {
        return false;
    }

    const double x = d.dot(t1);
    const double y = d.dot(t2);

    theta = std::atan2(y, x);
    if (theta < 0.0) {
        theta += 2.0 * M_PI;
    }

    return true;
}

// Returns true if two angular values are effectively the same direction.
bool anglesCoincident(double a, double b, double eps) {
    double diff = std::abs(a - b);
    diff = std::min(diff, 2.0 * M_PI - diff);
    return diff <= eps;
}

// Given two unit vectors, compute the rotation from one to the other
// Rotation formulation taken from https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
Eigen::Matrix3d computeRotation(const Eigen::Vector3d& u, const Eigen::Vector3d& v) {
    const double eps = 1e-8;

    if (u.norm() <= eps || v.norm() <= eps) {
        return Eigen::Matrix3d::Identity();
    }

    Eigen::Vector3d a = u.normalized();
    Eigen::Vector3d b = v.normalized();

    double cosUV = std::clamp(a.dot(b), -1.0, 1.0);

    if (cosUV >= 1.0 - eps) {
        return Eigen::Matrix3d::Identity();
    }

    if (cosUV <= -1.0 + eps) {
        Eigen::Vector3d axis = a.unitOrthogonal();
        return Eigen::AngleAxisd(M_PI, axis).toRotationMatrix();
    }

    Eigen::Vector3d axis = a.cross(b).normalized();
    double theta = std::acos(cosUV);

    return Eigen::AngleAxisd(theta, axis).toRotationMatrix();
}

// Overload
Eigen::Matrix3d computeRotation(const Eigen::Vector3d& axis, const double& theta) {
    const double eps = 1e-8;

    if (axis.norm() <= eps) {
        return Eigen::Matrix3d::Identity();
    }

    return Eigen::AngleAxisd(theta, axis.normalized()).toRotationMatrix();
}

// Find basis vectors for a planar region (ex. tangent plane)
// Build plane basis given only n, and unitialized t1, t2
void buildPlaneBasis(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2) {
    const double eps = 1e-12;
    if (n.squaredNorm() <= eps) {
        t1 = Eigen::Vector3d::UnitX();
        t2 = Eigen::Vector3d::UnitY();
        return;
    }

    Eigen::Vector3d n_norm = n.normalized();
    if (std::abs(n(0)) < 0.9)
        t1 = n_norm.cross(Eigen::Vector3d::UnitX()).normalized();
    else
        t1 = n_norm.cross(Eigen::Vector3d::UnitY()).normalized();
    t2 = n_norm.cross(t1); // already unit
    return;
}

// Given a point on a plane basis and the plane basis, convert to 2D planar point
Eigen::Vector2d convertTo2D(const Eigen::Vector3d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
    Eigen::Vector3d vec = p - origin;
    return Eigen::Vector2d(vec.dot(t1), vec.dot(t2));
}

// Given a 2D planar point and the plane basis, revert to its 3D counterapart
Eigen::Vector3d revertTo3D(const Eigen::Vector2d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
    return origin + p(0) * t1 + p(1) * t2;
}

// Projects a vector onto the tangent plane of a normal vector.
// If the resulting projection is near-degenerate, then return a random unit tangent vector.
double projectVectorOntoTangentPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& vec, Eigen::Vector3d& proj, double scale) {
    const double eps = 1e-12;
    // Ensure normal is unit
    Eigen::Vector3d n = normal.normalized();

    // Project onto tangent plane (i.e., subtract out projection onto normal vector)
    proj = vec - vec.dot(n) * n;
    double len = proj.norm();

    // Under degeneracy, sample a random unit vector in the tangent plane using some
    // constructed tangent plane basis
    if (len < eps) {
        // Build orthonormal tangent basis
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        buildPlaneBasis(normal, t1, t2);

        // Pick vector by sampling a random angle in [0, 2pi)
        static std::mt19937 gen(std::random_device{}());
        static std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
        // Construct proj using basis vectors
        double theta = angle_dist(gen);
        proj = std::cos(theta) * t1 + std::sin(theta) * t2;
        return -1;
    }
    // Normalize, scale, and return length
    proj /= len;
    proj *= scale;
    return len;
}

// Projects a point onto the tangent plane of a normal given a center 
Eigen::Vector3d projectPointOntoPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& center, const Eigen::Vector3d& p) {
    const double eps = 1e-12;
    if (normal.squaredNorm() <= eps) {
        return p;
    }
    Eigen::Vector3d n = normal.normalized();
    return p - (p - center).dot(n) * n;
}

// Check if a 2D point is in a 2D polygon
// To do this, we do raycasting to the segment
bool pointInPolygon2D(const Eigen::Vector2d& p, const std::vector<Eigen::Vector2d>& poly) {
    bool inside = false;
    int n = poly.size();

    for (int v0 = 0; v0 < n; v0++) {
        int v1 = (v0 + n - 1) % n;
        const Eigen::Vector2d& p0 = poly[v0];
        const Eigen::Vector2d& p1 = poly[v1];
        // If 
        bool intersect = ((p0(1) > p(1)) != (p1(1) > p(1))) &&
                         (p(0) < (p1(0) - p0(0)) * (p(1) - p0(1)) / (p1(1) - p0(1)) + p0(0));

        if (intersect) {
            inside = !inside;
        }
    }
    return inside;
}

// Get the closest point on a segment in 2D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
Eigen::Vector2d closestPointOnSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, bool clip) {
    Eigen::Vector2d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16) {
        return v0;
    }
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return v0 + t * vec;
}

double cross2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

// TODO: Fix this function so that it properly sets u
bool raycastToSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& direc, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1,
                        double& t, double& u) {
    double eps = 1e-12;
    const Eigen::Vector2d seg = v1 - v0;
    Eigen::Vector2d direc_norm = direc.normalized();
    const double denom = cross2D(direc_norm, seg);

    // Parallel or nearly parallel
    if (std::abs(denom) < eps) {
        return false;
    }

    const Eigen::Vector2d rhs = v0 - p;

    t = cross2D(rhs, seg) / denom;
    u = cross2D(rhs, direc_norm) / denom;

    // Ray constraint and segment constraint
    if (t < -eps) {
        return false;
    }
    if (u < -eps || u > 1.0 + eps) {
        return false;
    }

    // Clamp tiny numerical drift
    if (t < 0.0) {
        t = 0.0;
    }
    if (u < 0.0) {
        u = 0.0;
    } else if (u > 1.0) {
        u = 1.0;
    }
    return true;
}

// Get the closest point on a segment in 3D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
Eigen::Vector3d closestPointOnSegment3D(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, bool clip) {
    Eigen::Vector3d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16)
        return v0;
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return v0 + t * vec;
}

// MEAN VALUE COORDINATES
// Compute angle between any two 2D vectors given the four endpoints
// Vectors are computed as p1 - p0, p3 - p2
double vectorAngle(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3) {
    Eigen::Vector2d v0 = p1 - p0;
    Eigen::Vector2d v1 = p3 - p2;
    return std::atan2(v0(0)*v1(1) - v0(1)*v1(0), v0(0)*v1(0) + v0(1)*v1(1));
}

// Compute the sign of a double value
double computeSign(const double& value) {
    if (value > 0.0) {
        return 1.0;
    } else if (value < 0.0) {
        return -1.0;
    }
    return 0.0;
}

void meanValueCoordinates(const Eigen::Vector2d& target, const std::vector<Eigen::Vector2d>& cage, Eigen::VectorXd& weights) {
    weights.setZero();

    double W = 0.0;
    std::vector<double> beta(cage.size());
    std::vector<double> gamma(cage.size());
    std::vector<double> s(cage.size());
    std::vector<double> r(cage.size());
    for (int v = 0; v < cage.size(); v++) {
        beta[v] = vectorAngle(cage[v], cage[(v+1)%cage.size()], target, cage[v]);
        gamma[v] = vectorAngle(cage[(v+1)%cage.size()], cage[v], cage[(v+1)%cage.size()], target);
        s[v] = beta[v] + gamma[v];
        r[v] = (cage[v] - target).norm();
    }
    std::vector<double> w(cage.size());
    for (int v = 0; v < cage.size(); v++) {
        int v_m1 = (v + cage.size() - 1)%cage.size();
        double alpha_m1p1 = vectorAngle(target, cage[v_m1], target, cage[(v+1)%cage.size()]);
        double s_m1p1 = M_PI * (computeSign(s[v_m1]) + computeSign(s[v])) - s[v_m1] - s[v];
        if (computeSign(alpha_m1p1) != computeSign(s_m1p1)) {    // NOTE: Potential numerical issue here
            alpha_m1p1 *= -1.0;
        }
        w[v] = r[v_m1] * std::sin(alpha_m1p1/2.0);
        for (int u = 0; u < cage.size(); u++) {
            if ((u != v_m1) && (u != v)) {
                w[v] *= r[u] * std::sin(std::abs(s[v])/2.0);
            }
        }
        W += w[v];
    }
    // Normalize to get final weights
    for (int v = 0; v < cage.size(); v++) {
        weights[v] = w[v]/W;
    }

    return;
}

} // namespace Utils