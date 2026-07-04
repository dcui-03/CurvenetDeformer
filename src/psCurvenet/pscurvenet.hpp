// pscurvenet.hpp
#pragma once

#include "pscurvenet_types.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace psCurvenet {

// Editable curvenet class with various operations for test editing
// NOTE: This curvenet representation is NOT associated with profilemover class and is for Polyscope testing purposes only
class pscurvenet {
    public:
        // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
        // NOTE: Constructor assumes you already have no duplicates in your inputs
        pscurvenet(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents, std::vector<std::array<int, 4>> Splines);
        // Empty constructor
        pscurvenet();

        // --------- UPDATE CURVENET -----------
        // Hard reset everything
        void resetCurvenet();
        // Update control position
        void updateControlPos(int c, Eigen::Vector3d new_pos);
        // Update tangent position
        void updateTangentPos(int t, Eigen::Vector3d new_pos);
        // Add a control and return its index
        int addControl(Eigen::Vector3d pos);
        // Add a spline and return its index
        int addSpline(int c0, Eigen::Vector3d t0_pos, Eigen::Vector3d t1_pos, int c1);
        // Remove a control. Requires some major reordering...
        int removeControl(int c);
        // Remove a spline. Also requires some reordering...
        int removeSpline(int s);

        // --------- GETTERS + POLYSCOPE CONVERSION -----------
        // Control positions but returned as an Eigen::MatrixXd
        Eigen::MatrixXd cPosAsMatrix();
        // Tangent positions but returned as an Eigen::MatrixXd
        Eigen::MatrixXd tPosAsMatrix();
        // Curve network as a chain of discrete splines
        void cnAsCurveNetwork(Eigen::MatrixXd& verts, std::vector<std::pair<int, int>> connectivity);
        // Curve network as a chain of discrete splines with extra segments for the tangents
        void cnAsCurveNetworkWithTans(Eigen::MatrixXd& verts, std::vector<std::pair<int, int>> connectivity);


        // --------- SAMPLING -----------
        // Sample a bezier curve at time t
        Eigen::Vector3d tSampleBezier(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const;
        Eigen::Vector3d tSampleBezier(int s, double t) const;
        // NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
        // Returns n_samples points on the curve, including the endpoints
        std::vector<Eigen::Vector3d> sampleBezierNaive(int s, int n_samples = 50) const;
    protected:
        // No class inheritance
    private:
        // --------- Internal Deletion -----------
        // To remove a control, first collect all its adjacent halfedges and splines
        // Then remove the control from the map and change the indices of the later controls in the map
        // Then do the same for the halfedges and splines
        // No this won't work..... What to do... disallow?
        int deleteControl(int c);
        int deleteHalfedges(int he);    // Delete a halfedge AND its twin
        int deleteSpline(int s);

        // --------- ITERATORS -----------
        // Get the adjacent splines to a control vertex
        std::vector<int> ctrlAdjSplines(int c) const;

        // Returns vertices adjacent to a spline
        std::vector<int> splineAdjCtrls(int s) const;

        // Store attributes as maps, where the key is the index
        std::map<int, Control> C;
        std::map<int, HalfEdge> HE;
        std::map<int, Spline> S;
        // Store next valid key
        int nextC = 0;
        int nextHE = 0;
        int nextS = 0;
};

}   // namespace psCurvenet