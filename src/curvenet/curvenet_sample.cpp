#include "curvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <cmath>


namespace Curvenet {
    // Sample a bezier curve at time t
    Eigen::Vector3d curvenet::tSampleBezier(Eigen::Vector3d c0, Eigen::Vector3d c1, Eigen::Vector3d c2, Eigen::Vector3d c3, double t) {
        Eigen::Vector3d sample = std::pow(1 - t, 3) * c0 +
                        3 * std::pow(1 - t, 2) * t * c1 +
                        3 * (1 - t) * std::pow(t, 2) * c2 +
                        std::pow(t, 3) * c3;
        return sample;
    }
    Eigen::Vector3d curvenet::tSampleBezier(int s, double t) {
        Eigen::Vector3d c0, c1, c2, c3;
        int he = S[s].he;
        c0 = C[HE[he].origin].pos;
        c1 = HE[he].tan;
        c3 = HE[HE[he].twin].tan;
        c3 = C[HE[HE[he].twin].origin].pos;
        return tSampleBezier(c0, c1, c2, c3, t);
    }

    // NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
    // Returns n_samples+1 points on the curve, including the endpoints
    std::vector<Eigen::Vector3d> curvenet::sampleBezierNaive(int s, int n_samples = 50) {
        std::vector<Eigen::Vector3d> samples(n_samples);
        int he = S[s].he;
        double h = 1.0/(n_samples - 1);
        double t = 0.0;

        samples[0] = C[HE[he].origin].pos;  // Start
        for (int i = 1; i < n_samples; i++) {
            t += h;
            t = std::min(1.0, t);
            samples[i] = tSampleBezier(s, t);
        }
        samples[n_samples-1] = C[HE[HE[he].twin].origin].pos;   // End

        return samples;
    }

    // Estimate the arclength
    double curvenet::arclenEst(std::vector<Eigen::Vector3d> samples) {
        double length = 0.0;
        for (int i = 0; i < samples.size() - 1; i++) {
            length += (samples[i+1] - samples[i]).norm();
        }
        return length;
    }
    double curvenet::arclenEst(int s, int n_samples = 50) {
        std::vector<Eigen::Vector3d> samples = sampleBezierNaive(s, n_samples);
        return arclenEst(samples);
    }

    // Uniformly sample based on arclength estimator
    std::vector<Eigen::Vector3d> curvenet::unifSample(int s, int n_samples) {
        // 1. Estimate the length of the bezier curve
        // TODO: Pick a number of samples for the regular sampling
        std::vector<Eigen::Vector3d> regSamples = sampleBezierNaive(s);
        double arclength = arclenEst(regSamples);

        // 2. Take discrete samples off of the regular sampling
        // Get the n desired samples
        double h = 1.0/(n_samples - 1);
        double t = 0.0;
        std::vector<Eigen::Vector3d> unifSamples(n_samples);

        // track the index of the next vertex
        int next_regV = 1;
        double next_regDist = (regSamples[next_regV] - regSamples[next_regV - 1]).norm();
        for (int i = 1; i < n_samples; i++) {
            // track where we are in our walk
            double h_local = h;
            while (h_local > 0.0) { // Walk until we've walked the right distance
                h_local -= next_regDist;
                if (next_regDist >= h_local) {  // Walks to the next edge
                    next_regV++;
                    if (next_regV == regSamples.size() - 1) {
                        break;
                    }
                    next_regDist = (regSamples[next_regV] - regSamples[next_regV - 1]).norm();
                } else {
                    next_regDist -= h_local;
                }
            }
            unifSamples[i] = regSamples[next_regV - 1] -  next_regDist * (regSamples[next_regV-1] - regSamples[next_regV]).normalized();
        }
        
        return unifSamples;
    }

}   // namespace Curvenet