#include "spline.hpp"

#include "control.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <cmath>


namespace Curvenet {

    spline::spline(std::vector<control>* cList, std::vector<tangent>* tList,
                   std::array<int, 2> controlIdx, std::array<int, 2> tangentIdx):
        controlList(cList), tangentList(tList),
        c0(controlIdx[0]), c1(controlIdx[1]), t0(tangentIdx[0]), t1(tangentIdx[1]) {
    }

    // Query the spline at a point t
    Eigen::Vector3d spline::computeSample(double t) {
        Eigen::Vector3d p0 = (controlList->at(c0)).getPos();
        Eigen::Vector3d p1 = (tangentList->at(t0)).getPos();
        Eigen::Vector3d p2 = (tangentList->at(t1)).getPos();
        Eigen::Vector3d p3 = (controlList->at(c1)).getPos();
        Eigen::Vector3d result = std::pow(1 - t, 3) * p0 +
                        3 * std::pow(1 - t, 2) * t * p1 +
                        3 * (1 - t) * std::pow(t, 2) * p2 +
                        std::pow(t, 3) * p3;
        return result;
    }

    // Compute a rule for how many samples to take as a function of user input, mean edge length, spline length
    // NOTE: You must compute edge length before 
    int spline::computeNumSamples(int alpha, double meanE, double arclength) {
        return std::max(2, static_cast<int>(alpha * (arclength)/meanE));
    }

    // Sample n vertices on the spline via arclength estimation
    // NOTE: Does NOT include endpoints
    std::vector<Eigen::Vector3d> spline::nSamplesByArclength(int alpha, double m_meanE, bool include_ends) {
        std::vector<Eigen::Vector3d> arclengthSamples;
        double arclength = UniformSampling(arclengthSamples, true, true);
        int n = computeNumSamples(alpha, m_meanE, arclength);

        // Accumulate lengths for easy traversal
        double total = 0.0;
        std::vector<double> arclengthPercents(arclengthSamples.size());
        arclengthPercents[0] = 0.0;
        arclengthPercents[arclengthSamples.size() - 1] = 1.0;
        for (int e = 1; e < arclengthSamples.size() - 1; e++) {
            total += (arclengthSamples[e+1] - arclengthSamples[e-1]).norm();
            arclengthPercents[e] = total / arclength;
        }
        
        // Get the n desired samples
        double h = 1.0/(n + 1);
        double t = 0.0;
        std::vector<Eigen::Vector3d> nSamples(n);
        for (int s = 0; s < n; s++) {
            t += h;
            // Check which "bracket" we're in
            for (int e = 0; e < arclengthSamples.size() - 1; e++) {
                if ((arclengthPercents[e] <= t) && (arclengthPercents[e+1] > t)) {
                    double t_sub;
                    if (t - arclengthPercents[e] < 1e-5) {
                        t_sub = 0.0;
                    } else {
                        t_sub = (t - arclengthPercents[e])/(arclengthPercents[e+1] - arclengthPercents[e]);
                    }
                    nSamples[s] = arclengthSamples[e] + t_sub * (arclengthSamples[e+1] - arclengthSamples[e]);
                    break;
                }
            }
        }
        if (include_ends) {
            nSamples.insert(nSamples.begin(), controlList->at(c0).getPos());
            nSamples.push_back((controlList->at(c1)).getPos());
        }
        return nSamples;
    }

    // Sample n points uniformly to get an estimated arclength
    // Adaptively sample based on some rules
    double spline::UniformSampling(std::vector<Eigen::Vector3d>& samplePoints, bool include_ends, bool adaptiveSampling) {
        samplePoints.clear();
        int numSamples = 40;
        if (adaptiveSampling) {
            numSamples = std::max(2, adaptiveSamplingRule());
        }
        double h = 1.0/(numSamples + 1);
        double t = 0.0;
        // Add the start
        if (include_ends) {
            samplePoints.push_back((controlList->at(c0)).getPos());
        }
        for (int i = 0; i < numSamples; i++) {
            t += h;
            t = std::min(1.0, t);
            samplePoints.push_back(computeSample(t));
        }
        // Add the end
        if (include_ends) {
            samplePoints.push_back((controlList->at(c1)).getPos());
        }

        // compute arclength
        double length = 0.0;
        if (include_ends) {
            for (int e = 1; e < samplePoints.size(); e++) {
                length += (samplePoints[e] - samplePoints[e-1]).norm();
            }
        } else {
            length += (samplePoints[0] - (controlList->at(c0)).getPos()).norm();
            length += (samplePoints[numSamples - 1] - (controlList->at(c1)).getPos()).norm();
            for (int e = 1; e < numSamples; e++) {
                length += (samplePoints[e] - samplePoints[e-1]).norm();
            }
        }
        return length;
    }

    // Number of samples to take for arclength estimation, using the control cage as a metric
    // TODO: Not sure what to do here... current set up is not arclength adaptive
    int spline::adaptiveSamplingRule(double meanE) {
        Eigen::Vector3d p0 = (controlList->at(c0)).getPos();
        Eigen::Vector3d p1 = (tangentList->at(t0)).getPos();
        Eigen::Vector3d p2 = (tangentList->at(t1)).getPos();
        Eigen::Vector3d p3 = (controlList->at(c1)).getPos();
        double v0 = (p1 - p0).norm();
        double v1 = (p2 - p1).norm();
        double v2 = (p3 - p2).norm();
        double enddiff = (p0 - p3).norm();
        double cagelen = v0 + v1 + v2;
        /*
        double angleWeight = std::min((p1 - p0).dot(p2 - p3), 0.0);
        if (angleWeight < 0.0) {
            angleWeight = 1.0 - (angleWeight);
        } else {
            angleWeight = 1.0;
        }*/
        return static_cast<int>((cagelen/enddiff) * 32);
    }


}   // namespace Curvenet