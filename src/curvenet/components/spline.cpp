#include "spline.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Curvenet {

spline::spline(std::vector<int> splineP) {
    if (splineP.size() != 4) {
        throw std::invalid_argument("Spline needs 4 indices [start,t1,t2,end]");
    }
    c1 = splineP[0];
    t1 = splineP[1];
    t2 = splineP[2];
    c2 = splineP[3];
}

void spline::computeSpline(const std::vector<control>& controls, const std::vector<tangent>& tangents) {
    if (c1 < 0 || c2 < 0 || t1 < 0 || t2 < 0) {
        throw std::runtime_error("Spline indices are not initialized");
    }

    const std::size_t c1Idx = static_cast<std::size_t>(c1);
    const std::size_t c2Idx = static_cast<std::size_t>(c2);
    const std::size_t t1Idx = static_cast<std::size_t>(t1);
    const std::size_t t2Idx = static_cast<std::size_t>(t2);

    if (c1Idx >= controls.size() || c2Idx >= controls.size() || t1Idx >= tangents.size() ||
        t2Idx >= tangents.size()) {
        throw std::out_of_range("Spline references out-of-range control/tangent index");
    }

    p0 = controls[c1Idx].getPosition();
    p1 = controls[c2Idx].getPosition();
    h0 = tangents[t1Idx].getPosition();
    h1 = tangents[t2Idx].getPosition();
    arclength = computeArclength();
}

int spline::computeNumSamples(int alpha, double meanE) const {
    if (alpha <= 0) {
        return 2;
    }
    if (meanE <= std::numeric_limits<double>::epsilon()) {
        return std::max(2, alpha);
    }
    const int byLength = static_cast<int>(std::ceil(arclength / meanE));
    return std::max(2, alpha * std::max(1, byLength));
}

std::vector<Eigen::Vector3d> spline::sampleParameterization(int n) const {
    if (n < 2) {
        n = 2;
    }

    std::vector<Eigen::Vector3d> out;
    out.reserve(static_cast<std::size_t>(n));

    const bool hasArc = !arcLength.empty();
    for (int i = 0; i < n; ++i) {
        const double u = static_cast<double>(i) / static_cast<double>(n - 1);
        if (hasArc && arclength > std::numeric_limits<double>::epsilon()) {
            out.push_back(tSampleParameterization(parameterAtArclength(u * arclength)));
        } else {
            out.push_back(tSampleParameterization(u));
        }
    }
    return out;
}

std::vector<Eigen::Vector3d> spline::sampleUniformParameterization(int n) const {
    if (n < 2) {
        n = 2;
    }

    std::vector<Eigen::Vector3d> out;
    out.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(n - 1);
        out.push_back(tSampleParameterization(t));
    }
    return out;
}

std::vector<Eigen::Vector3d> spline::nSamplesByArclength(int alpha, double meanE, bool include_ends) const {
    std::vector<Eigen::Vector3d> samples = sampleParameterization(computeNumSamples(alpha, meanE));
    if (include_ends || samples.size() <= 2) {
        return samples;
    }
    return std::vector<Eigen::Vector3d>(samples.begin() + 1, samples.end() - 1);
}

double spline::UniformSampling(std::vector<Eigen::Vector3d>& samplePoints,
                               bool include_ends,
                               bool adaptiveSampling) const {
    const int interiorSamples = adaptiveSampling ? adaptiveSamplingRule() : 40;
    const std::vector<Eigen::Vector3d> full = sampleUniformParameterization(std::max(2, interiorSamples + 2));

    if (include_ends || full.size() <= 2) {
        samplePoints = full;
    } else {
        samplePoints.assign(full.begin() + 1, full.end() - 1);
    }

    double length = 0.0;
    for (std::size_t i = 1; i < full.size(); ++i) {
        length += (full[i] - full[i - 1]).norm();
    }
    return length;
}

double spline::controlPolylineLength() const {
    return (p0 - h0).norm() + (h0 - h1).norm() + (h1 - p1).norm();
}

Eigen::Vector3d spline::outgoingTangentAtControl(int controlIdx) const {
    if (controlIdx == c1) {
        return h0 - p0;
    }
    if (controlIdx == c2) {
        return h1 - p1;
    }
    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d spline::tSampleParameterization(double t) const {
    t = clamp01(t);
    const double u = 1.0 - t;
    return p0 * (u * u * u) + h0 * (3.0 * u * u * t) + h1 * (3.0 * u * t * t) + p1 * (t * t * t);
}

double spline::computeArclength() {
    constexpr std::size_t kSamples = 256;

    arcT.assign(kSamples, 0.0);
    arcLength.assign(kSamples, 0.0);

    double cumulative = 0.0;
    Eigen::Vector3d prev = tSampleParameterization(0.0);

    for (std::size_t i = 1; i < kSamples; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(kSamples - 1);
        const Eigen::Vector3d curr = tSampleParameterization(t);
        cumulative += (curr - prev).norm();
        arcT[i] = t;
        arcLength[i] = cumulative;
        prev = curr;
    }

    return cumulative;
}

double spline::parameterAtArclength(double s) const {
    if (arcLength.empty() || arcT.empty() || arclength <= std::numeric_limits<double>::epsilon()) {
        return 0.0;
    }
    s = std::clamp(s, 0.0, arclength);

    auto it = std::lower_bound(arcLength.begin(), arcLength.end(), s);
    if (it == arcLength.begin()) {
        return 0.0;
    }
    if (it == arcLength.end()) {
        return 1.0;
    }

    const std::size_t i1 = static_cast<std::size_t>(it - arcLength.begin());
    const std::size_t i0 = i1 - 1;
    const double s0 = arcLength[i0];
    const double s1 = arcLength[i1];
    const double t0 = arcT[i0];
    const double t1v = arcT[i1];

    if (std::abs(s1 - s0) <= std::numeric_limits<double>::epsilon()) {
        return t0;
    }

    const double w = (s - s0) / (s1 - s0);
    return t0 + w * (t1v - t0);
}

int spline::adaptiveSamplingRule(double meanE) const {
    const double enddiff = (p0 - p1).norm();
    const double cageLen = controlPolylineLength();

    int rule = 16;
    if (enddiff > std::numeric_limits<double>::epsilon()) {
        rule = std::max(rule, static_cast<int>(std::ceil((cageLen / enddiff) * 32.0)));
    }
    if (meanE > std::numeric_limits<double>::epsilon()) {
        rule = std::max(rule, static_cast<int>(std::ceil(cageLen / meanE)));
    }
    return std::max(2, rule);
}

double spline::clamp01(double t) { return std::clamp(t, 0.0, 1.0); }

}  // namespace Curvenet
