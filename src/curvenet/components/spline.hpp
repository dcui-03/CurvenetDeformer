// spline.hpp
#pragma once

#include "control.hpp"

#include <Eigen/Core>

#include <vector>

namespace Curvenet {

class spline {
public:
    spline() = default;
    explicit spline(std::vector<int> splineP);

    void computeSpline(const std::vector<control>& controls, const std::vector<tangent>& tangents);

    int computeNumSamples(int alpha, double meanE) const;
    std::vector<Eigen::Vector3d> sampleParameterization(int n) const;
    std::vector<Eigen::Vector3d> sampleUniformParameterization(int n) const;

    // Legacy editing helpers preserved for the older UI path.
    std::vector<Eigen::Vector3d> nSamplesByArclength(int alpha, double meanE, bool include_ends = true) const;
    double UniformSampling(std::vector<Eigen::Vector3d>& samplePoints,
                           bool include_ends = true,
                           bool adaptiveSampling = true) const;

    double controlPolylineLength() const;
    Eigen::Vector3d outgoingTangentAtControl(int controlIdx) const;

    int startControl() const { return c1; }
    int endControl() const { return c2; }
    int startTangent() const { return t1; }
    int endTangent() const { return t2; }

private:
    Eigen::Vector3d tSampleParameterization(double t) const;
    double computeArclength();
    double parameterAtArclength(double s) const;
    int adaptiveSamplingRule(double meanE = -1.0) const;
    static double clamp01(double t);

    int c1 = -1;
    int c2 = -1;
    int t1 = -1;
    int t2 = -1;
    double arclength = 0.0;

    Eigen::Vector3d p0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d p1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d h0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d h1 = Eigen::Vector3d::Zero();

    std::vector<double> arcT;
    std::vector<double> arcLength;
};

}  // namespace Curvenet
