// curvenet.hpp
#pragma once

#include "curvenet_types.hpp"
#include <Eigen/Core>
#include <array>

namespace Curvenet {

// Curvenet class consisting of cubic bezier splines
// NOTE: This class is NOT meant to be an editable curve network. It is purely static and has the purpose of supporting a profilemover object.
class curvenet {
    public:
        // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
        // NOTE: Constructor assumes you already have no duplicates in your inputs
        curvenet(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents, std::vector<std::array<int, 4>> Splines);
        // Empty constructor
        curvenet();

        // --------- SAMPLING -----------
        // Sample a bezier curve at time t
        Eigen::Vector3d tSampleBezier(Eigen::Vector3d c0, Eigen::Vector3d c1, Eigen::Vector3d c2, Eigen::Vector3d c3, double t);
        Eigen::Vector3d tSampleBezier(int s, double t);
        // NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
        // Returns n_samples+1 points on the curve, including the endpoints
        std::vector<Eigen::Vector3d> sampleBezierNaive(int s, int n_samples = 50);
        // Estimate the arclength
        double arclenEst(int s, int n_samples = 50);
        double arclenEst(std::vector<Eigen::Vector3d>);
        // Uniformly sample based on arclength estimator
        // Takes a user parameter alpha which helps control sampling
        std::vector<Eigen::Vector3d> unifSample(int s, int n_samples = 50);

        // --------- EDITING -----------
        // Exposed position edits
        int editControlPos(int c, const Eigen::Vector3d pos);
        // Exposed normal augmentation
        // NOTE: Normals should only applied to vertices once, by projection onto the mesh
        int editControlN(int c, Eigen::Vector3d normal);

        // --------- OTHER -----------
        int sortAdjHEAll();
        int assignCtrlTypeAll(int c);
    protected:
        // No class inheritance
    private:
        // --------- INITIALIZATION -----------
        // Create new control
        int addControl(Eigen::Vector3d pos);

        // Add a spline given the start, end, and two tangent endpoints
        int addSpline(int start, int end, Eigen::Vector3d t0, Eigen::Vector3d t1);

        // --------- ITERATORS -----------
        // Get the adjacent tangent vectors to a control vertex
        std::vector<Eigen::Vector3d> ctrlAdjTans(int c);

        // Returns a list of the splines adjacent to a control 
        // If no self loops, then these are CCW
        std::vector<int> ctrlAdjSplines(int c);

        // Returns vertices adjacent to a spline
        std::vector<int> curvenet::splineAdjCtrls(int s);

        //  --------- OTHER -----------
        // Sort the halfedges of a control to be CCW
        // Should only be performed AFTER assigning normals to all verts
        int sortAdjHE(int c);
        // Compute what kind of vertex each control is using the valence of splines
        int assignCtrlType(int c);

        // Store attributes as lists
        std::vector<Control> C;
        std::vector<HalfEdge> HE;
        std::vector<CubicSpline> S;
        
};

}   // namespace Curvenet