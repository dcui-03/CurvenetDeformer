// spline.hpp
#pragma once

#include "control.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>

// A class for cubic Bézier splines

namespace Curvenet {

class spline {
    public:
        // Constructor takes four indices that the curvenet will make use of [start, tangent 1, tangent 2, end]
        spline(std::vector<control>* cList, std::vector<tangent>* tList,
                   std::array<int, 2> controlIdx, std::array<int, 2> tangentIdx);

        // Sample n points via arclength using adaptive sampling
        // NOTE: Precompute m_meanE as the user input / mean edge length in profile mover class
        std::vector<Eigen::Vector3d> nSamplesByArclength(int alpha, double m_meanE, bool include_ends = true);

        // Sample n points uniformly to get an estimated arclength
        double UniformSampling(std::vector<Eigen::Vector3d>& samplePoints, bool include_ends = true, bool adaptiveSampling = true);
        
        // Store indices of the associated controls and tangents in their lists
        int c0;
        int c1;
        int t0;
        int t1;
    protected:
        // No class inheritance
    private:
        // Compute a rule for how many samples to take as a function of user input, mean edge length, spline length
        int computeNumSamples(int alpha, double meanE, double arclength);

        // Number of samples to take for arclength estimation, using the control cage as a metric
        int adaptiveSamplingRule(double meanE = -1.0);

        // Compute own arclength
        // TODO: remove?
        double computeArclength(int n_samples);


        // Compute a sample on the spline at instance t in [0, 1]
        Eigen::Vector3d computeSample(double t);

        // Pointer to control list
        std::vector<control> *controlList;
        // Pointer to tangent list
        std::vector<tangent> *tangentList;
        // double arclength;
          
};

}   // namespace Curvenet