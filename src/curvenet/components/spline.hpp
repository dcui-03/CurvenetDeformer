// spline.hpp
#pragma once

#include "control.hpp"
#include <Eigen/Core>

#include <cstddef>
#include <vector>


namespace Curvenet {

class spline {
public:

    // Constructor takes four indices that the curvenet will make use of [start, tangent 1, tangent 2, end]
    spline(std::vector<int> splineP);

    // Compute the spline --> store local endpoints/tangents from parent arrays
    void computeSpline(const std::vector<control>& controls, const std::vector<tangent>& tangents);

    // Compute a rule for how many samples to take as a function of user input, mean edge length, spline length
    int computeNumSamples(int alpha, double meanE) const;

    // Sample n points via arclength.
    std::vector<Eigen::Vector3d> sampleParameterization(int n) const;

    // Uniform in parameter t in [0,1].
    std::vector<Eigen::Vector3d> sampleUniformParameterization(int n) const;

    // Length of the control polygon p0-h0-h1-p1.
    double controlPolylineLength() const;

    // Tangent direction pointing from a control toward the spline interior.
    Eigen::Vector3d outgoingTangentAtControl(int controlIdx) const;

    int startControl() const { return c1; }
    int endControl() const { return c2; }

protected:
    // No class inheritance
private:
    // Sample the point at t in range [0, 1] of the arclength parameterization
    Eigen::Vector3d tSampleParameterization(double t) const;

    // Compute own arclength
    double computeArclength();
    double parameterAtArclength(double s) const;
    static double clamp01(double t);

    // Store indices of the associated controls and tangents in their lists
    int c1 = -1;
    int c2 = -1;
    int t1 = -1;
    int t2 = -1;
    double arclength = 0.0;

    // Cached local spline geometry for sampling
    Eigen::Vector3d p0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d p1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d h0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d h1 = Eigen::Vector3d::Zero();

    std::vector<double> arcT;
    std::vector<double> arcLength;
};

}  // namespace Curvenet
