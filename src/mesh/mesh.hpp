// mesh.hpp
#pragma once

#include "dcurvenet/pdcurvenet.hpp"
#include "polyHE/polyHE.hpp"
#include <Eigen/Core>
#include <vector>


namespace Mesh {

// Template mesh class that augments our standard mesh setup with some quantities we need
class mesh {
    public:
        // Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
        // TODO: Add in the correct mesh data struct
        mesh(std::vector<Eigen::Vector3d>& V, std::vector<std::vector<int>>& F);

        // Project a vertex onto the mesh. If multiple, just picks the first one
        // Returns element type that was landed on (0 for vert, 1 for edge, 2 for face)
        // NOTE: Depending on the rule for nonplanar faces you may need to write your own version of this method
        int computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, int& elIdx, bool snap = true);

        // Compute the length of the diagonal of the bounding box.
        void computeBBoxDiag();

        // Given a point in space, computes the index of the nearest edge to that point on a given face
        double computeNearestFaceEdge(const int face, const Eigen::Vector3d& p, int& nearestIdx, Eigen::Vector3d& nearestPnt);
        
        // Getters
        Eigen::Vector3d& getVNormal(int vidx);
        Eigen::Vector3d& getENormal(int eidx);
        Eigen::Vector3d& getFNormal(int fidx);
        // Get pointer to he mesh
        polyHE::polyHE_t& getHEMesh();

        // Get mean edge length
        double getMeanE();
        // Get bbox diagonal
        double getBBoxDiag();

        int num_f;
        int num_v;

    protected:
        // Internal function to precompute height functions on both planar/nonplanar faces and their convexity
        virtual void computeHeightFuncsAndConvexity();
        // Internal function to precompute normals and areas on all mesh structures
        virtual void computeFNormalsAreas();
        // weight_fN weights by area
        virtual void computeENormals(bool weight_fN = false); // Note boundary edges get the adjacent face's normal.
        virtual void computeVNormalsAreas(bool weight_fN = true); // Average the nearby face normals using some weight. CAREFUL NOT TO GET DUPLICATE FACES with  cut-mesh
        
        // Computes mean edge length on the mesh
        void computeMeanE();

        // Mesh as a HE data structure
        polyHE::polyHE_t heMesh;
        // Neutral mesh state
        // Note: Do we need to store edge/face data or can we just let polyHE handle it?
        std::vector<Eigen::Vector3d>& V;
        std::vector<std::array<int, 2>> E;  // TODO: remove?
        std::vector<std::vector<int>>& F;
        std::vector<Eigen::VectorXd> H; // Height functions on faces as a VectorXd (let's us compute MVC result using a dot prod)
        std::vector<bool> Convex;   // Convexity of faces in the function (used for querying for Straightest Geodesic)

        // Local areas
        std::vector<double> vAreas;
        std::vector<double> fAreas;
        
        // Storage vectors by index for each normal
        std::vector<Eigen::Vector3d> vNormals;
        std::vector<Eigen::Vector3d> eNormals;
        std::vector<Eigen::Vector3d> fNormals;

        // Mean edge length on mesh
        double meanE;
        // Bounding box diagonal length
        // NOTE: For simplicity, we are using AABB
        double bboxDiag;
    private:
        // No unique structures
};

}   // namespace Mesh