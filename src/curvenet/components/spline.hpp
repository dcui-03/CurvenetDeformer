// spline.hpp
#pragma once

#include "control.hpp"
#include <Eigen/Core>
#include <vector>


namespace Curvenet {

class spline {
    public:
        // Constructor takes four indices that the curvenet will make use of [start, tangent 1, tangent 2, end]
        spline(std::vector<int> splineP);

        // Compute the spline --> Some function
        void computeSpline(std::vector<control>& controls, std::vector<tangents>& tangents);

        // Compute a rule for how many samples to take as a function of user input, mean edge length, spline length
        int computeNumSamples(int alpha, double meanE);

        // Sample n points via arclength.
        std::vector<Eigen::Vector3d> sampleParameterization(int n);
    protected:
        // No class inheritance
    private:
        // Sample the point at t in range [0, 1] of the arclength parameterization
        Eigen::Vector3d tSampleParameterization(double t);

        // Compute own arclength
        double computeArclength();

        // Store indices of the associated controls and tangents in their lists
        int c1;
        int c2;
        int t1;
        int t2;
        double arclength
        
        
};

}   // namespace CutMesh