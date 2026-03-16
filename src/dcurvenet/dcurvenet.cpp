#include "dcurvenet.hpp"
#include "utils/utils.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace DCurvenet {

namespace {

struct SideInterpolationResult {
    std::vector<Eigen::Vector3d> normals;
    std::vector<double> widths;
};

SideInterpolationResult interpolateSideAlongPolyline(const std::vector<Eigen::Vector3d>& tangents,
                                                     const std::vector<double>& lengths,
                                                     const Eigen::Vector3d& startNormal,
                                                     double startWidth,
                                                     bool blendToEnd,
                                                     const Eigen::Vector3d& endNormal,
                                                     double endWidth) {
    constexpr double kEps = 1e-12;
    constexpr double kWidthEps = 1e-8;

    const std::size_t k = tangents.size();
    SideInterpolationResult out;
    out.normals.assign(k, Eigen::Vector3d::Zero());
    out.widths.assign(k, std::max(kWidthEps, startWidth));
    if (k == 0) {
        return out;
    }

    std::vector<Eigen::Vector3d> omegaN(k, Eigen::Vector3d::Zero());
    Utils::projectVectorOntoTangentPlane(tangents[0], startNormal, omegaN[0]);
    for (std::size_t i = 1; i < k; ++i) {
        const Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(
            tangents[i - 1].normalized(),
            tangents[i].normalized());
        const Eigen::Vector3d rotated = q * omegaN[i - 1];
        Utils::projectVectorOntoTangentPlane(tangents[i], rotated, omegaN[i]);
    }

    double totalLen = 0.0;
    std::vector<double> prefixBefore(k, 0.0);
    for (std::size_t i = 0; i < k; ++i) {
        prefixBefore[i] = totalLen;
        totalLen += std::max(0.0, lengths[i]);
    }

    double theta = 0.0;
    if (blendToEnd) {
        Eigen::Vector3d endTarget;
        Utils::projectVectorOntoTangentPlane(tangents.back(), endNormal, endTarget);
        const Eigen::Vector3d tK = tangents.back().normalized();
        const Eigen::Vector3d nK = omegaN.back();
        const double num = nK.dot(endTarget.cross(tK));
        const double den = nK.dot(endTarget);
        theta = std::atan2(num, den);
    }

    for (std::size_t i = 0; i < k; ++i) {
        const Eigen::Vector3d ti = tangents[i].normalized();
        const double alpha = (totalLen > kEps) ? (prefixBefore[i] / totalLen) : 0.0;
        const Eigen::AngleAxisd torsion(alpha * theta, ti);
        const Eigen::Vector3d rotatedN = torsion * omegaN[i];
        Utils::projectVectorOntoTangentPlane(ti, rotatedN, out.normals[i]);

        if (blendToEnd) {
            out.widths[i] = std::max(kWidthEps, (1.0 - alpha) * startWidth + alpha * endWidth);
        } else {
            out.widths[i] = std::max(kWidthEps, startWidth);
        }
    }

    return out;
}

}  // namespace

dcurvenet::dcurvenet(const Curvenet::curvenet& CN,
                     double meanEdgeLength,
                     int samplesPerMeanEdge,
                     int uniformRefineSamples)
    : mean_edge_length(std::max(meanEdgeLength, std::numeric_limits<double>::epsilon())),
      samples_per_mean_edge(std::max(1, samplesPerMeanEdge)),
      uniform_refine_samples(std::max(8, uniformRefineSamples)) {
    initializeFromCurvenet(CN);
}

int dcurvenet::computeInteriorSampleCount(const Curvenet::spline& s,
                                          double meanEdgeLength,
                                          int samplesPerMeanEdge) {
    const double polylineLength = s.controlPolylineLength();
    const double scaled = static_cast<double>(std::max(1, samplesPerMeanEdge)) * (polylineLength / meanEdgeLength);
    return std::max(1, static_cast<int>(std::ceil(scaled)));
}

std::vector<Eigen::Vector3d> dcurvenet::resamplePolylineEvenArcLength(const std::vector<Eigen::Vector3d>& polyline,
                                                                       int outputCount) {
    if (polyline.empty()) {
        return {};
    }
    if (outputCount <= 1 || polyline.size() == 1) {
        return {polyline.front()};
    }

    std::vector<double> cumulative(polyline.size(), 0.0);
    for (std::size_t i = 1; i < polyline.size(); ++i) {
        cumulative[i] = cumulative[i - 1] + (polyline[i] - polyline[i - 1]).norm();
    }

    const double totalLength = cumulative.back();
    if (totalLength <= std::numeric_limits<double>::epsilon()) {
        std::vector<Eigen::Vector3d> collapsed(static_cast<std::size_t>(outputCount), polyline.front());
        collapsed.back() = polyline.back();
        return collapsed;
    }

    std::vector<Eigen::Vector3d> resampled;
    resampled.reserve(static_cast<std::size_t>(outputCount));

    for (int i = 0; i < outputCount; ++i) {
        const double target = (static_cast<double>(i) / static_cast<double>(outputCount - 1)) * totalLength;

        auto it = std::lower_bound(cumulative.begin(), cumulative.end(), target);
        if (it == cumulative.begin()) {
            resampled.push_back(polyline.front());
            continue;
        }
        if (it == cumulative.end()) {
            resampled.push_back(polyline.back());
            continue;
        }

        const std::size_t hi = static_cast<std::size_t>(it - cumulative.begin());
        const std::size_t lo = hi - 1;

        const double l0 = cumulative[lo];
        const double l1 = cumulative[hi];
        const double denom = l1 - l0;
        const double w = (denom <= std::numeric_limits<double>::epsilon()) ? 0.0 : (target - l0) / denom;

        resampled.push_back((1.0 - w) * polyline[lo] + w * polyline[hi]);
    }

    return resampled;
}

bool dcurvenet::outgoingAtControl(int controlDvertIdx,
                                  int dsplineIdx,
                                  int& segIdx,
                                  Eigen::Vector3d& tOut,
                                  double& segLen) const {
    segIdx = -1;
    tOut = Eigen::Vector3d::Zero();
    segLen = 0.0;

    if (dsplineIdx < 0 || static_cast<std::size_t>(dsplineIdx) >= dSplines.size()) {
        return false;
    }

    const dspline& ds = dSplines[static_cast<std::size_t>(dsplineIdx)];
    double sign = 1.0;

    if (ds.startDvert() == controlDvertIdx) {
        segIdx = ds.firstSegment();
        sign = 1.0;
    } else if (ds.endDvert() == controlDvertIdx) {
        segIdx = ds.lastSegment();
        sign = -1.0;
    } else {
        return false;
    }

    if (segIdx < 0 || static_cast<std::size_t>(segIdx) >= dSegments.size()) {
        return false;
    }

    const dsegment& seg = dSegments[static_cast<std::size_t>(segIdx)];
    segLen = seg.segmentLength();
    const Eigen::Vector3d tangent = sign * seg.direction();
    tOut = tangent.normalized();
    return true;
}

const SegmentData* dcurvenet::findControlSegmentData(int controlDvertIdx, int segmentIdx) const {
    if (controlDvertIdx < 0 || static_cast<std::size_t>(controlDvertIdx) >= dVerts.size()) {
        return nullptr;
    }
    const auto& data = dVerts[static_cast<std::size_t>(controlDvertIdx)].getSegmentData();
    for (const auto& d : data) {
        if (d.segment_idx == segmentIdx) {
            return &d;
        }
    }
    return nullptr;
}

void dcurvenet::computeControlRibbonData(const std::vector<Curvenet::control>& controls,
                                         const std::vector<int>& controlToDvert) {
    constexpr double kEps = 1e-12;
    constexpr double kWidthEps = 1e-8;

    for (std::size_t cIdx = 0; cIdx < controls.size(); ++cIdx) {
        const int dIdx = controlToDvert[cIdx];
        if (dIdx < 0 || static_cast<std::size_t>(dIdx) >= dVerts.size()) {
            continue;
        }

        dvert& dv = dVerts[static_cast<std::size_t>(dIdx)];
        const std::vector<int>& outgoing = dv.dsplineIdxs();
        const std::size_t n = outgoing.size();

        const Eigen::Vector3d surfaceN = controls[cIdx].getNormal().normalized();

        std::vector<SegmentData> segmentData;
        segmentData.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            SegmentData sd;
            sd.dspline_idx = outgoing[i];
            outgoingAtControl(dIdx, outgoing[i], sd.segment_idx, sd.t_out, sd.l);
            segmentData.push_back(sd);
        }
        dv.setSegmentData(segmentData);

        if (n == 0) {
            dv.setCornerNormals({});
            continue;
        }

        if (n <= 2) {
            dv.setCornerNormals(std::vector<Eigen::Vector3d>(n, surfaceN));

            for (std::size_t i = 0; i < n; ++i) {
                SegmentData& sd = segmentData[i];
                sd.n_plus = surfaceN;
                sd.n_minus = surfaceN;
                const double w = std::max(kWidthEps, sd.l);
                sd.w_plus = w;
                sd.w_minus = w;
            }
            dv.setSegmentData(segmentData);
            continue;
        }

        std::vector<Eigen::Vector3d> tangents(n, Eigen::Vector3d::Zero());
        for (std::size_t i = 0; i < n; ++i) {
            tangents[i] = segmentData[i].t_out;
        }

        std::vector<Eigen::Vector3d> corners(n, Eigen::Vector3d::Zero());
        std::vector<double> cNorm(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const Eigen::Vector3d& t0 = tangents[i];
            const Eigen::Vector3d& t1 = tangents[(i + 1) % n];
            corners[i] = t0.cross(t1);
            cNorm[i] = corners[i].norm();
        }

        std::vector<Eigen::Vector3d> m(n, Eigen::Vector3d::Zero());
        for (std::size_t i = 0; i < n; ++i) {
            Eigen::Vector3d cornerNormal = corners[i].normalized();
            if (cornerNormal.squaredNorm() <= kEps) {
                const Eigen::Vector3d blend = corners[(i + n - 1) % n] + corners[(i + 1) % n];
                cornerNormal = blend.normalized();
            }

            // Keep orientation consistent with the local surface normal.
            if (cornerNormal.dot(surfaceN) < 0.0) {
                cornerNormal *= -1.0;
            }
            m[i] = cornerNormal;
        }
        dv.setCornerNormals(m);

        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t prev = (i + n - 1) % n;
            const std::size_t next = (i + 1) % n;
            SegmentData& sd = segmentData[i];

            if (sd.segment_idx < 0 || static_cast<std::size_t>(sd.segment_idx) >= dSegments.size()) {
                continue;
            }

            const double li = segmentData[i].l;
            const double lPrev = segmentData[prev].l;
            const double lNext = segmentData[next].l;

            sd.n_plus = m[i];
            sd.n_minus = m[prev];
            sd.w_plus = std::max(kWidthEps, li + cNorm[i] * (lNext - li));
            sd.w_minus = std::max(kWidthEps, li + cNorm[prev] * (lPrev - li));
        }
        dv.setSegmentData(segmentData);
    }
}

void dcurvenet::interpolateRibbonAlongDSplines() {
    constexpr double kWidthEps = 1e-8;

    for (const dspline& ds : dSplines) {
        const auto& segIdxs = ds.segmentIdxs();
        const std::size_t k = segIdxs.size();
        if (k == 0) {
            continue;
        }

        std::vector<Eigen::Vector3d> t(k, Eigen::Vector3d::Zero());
        std::vector<double> l(k, 0.0);
        for (std::size_t i = 0; i < k; ++i) {
            const int segIdx = segIdxs[i];
            if (segIdx < 0 || static_cast<std::size_t>(segIdx) >= dSegments.size()) {
                continue;
            }
            t[i] = dSegments[static_cast<std::size_t>(segIdx)].direction();
            l[i] = dSegments[static_cast<std::size_t>(segIdx)].segmentLength();
        }

        const int startDv = ds.startDvert();
        const int endDv = ds.endDvert();
        const bool startIntersection =
            (startDv >= 0 && static_cast<std::size_t>(startDv) < dVerts.size() &&
             dVerts[static_cast<std::size_t>(startDv)].numOutgoing() >= 3);
        const bool endIntersection =
            (endDv >= 0 && static_cast<std::size_t>(endDv) < dVerts.size() &&
             dVerts[static_cast<std::size_t>(endDv)].numOutgoing() >= 3);

        const SegmentData* startSD = findControlSegmentData(startDv, segIdxs.front());
        const SegmentData* endSD = findControlSegmentData(endDv, segIdxs.back());
        if (startSD == nullptr || endSD == nullptr) {
            continue;
        }

        std::vector<Eigen::Vector3d> plusN(k, Eigen::Vector3d::Zero());
        std::vector<Eigen::Vector3d> minusN(k, Eigen::Vector3d::Zero());
        std::vector<double> plusW(k, std::max(kWidthEps, startSD->w_plus));
        std::vector<double> minusW(k, std::max(kWidthEps, startSD->w_minus));

        if (startIntersection || !endIntersection) {
            // Forward interpolation from start endpoint.
            const bool blendEnd = startIntersection && endIntersection;
            const SideInterpolationResult plus = interpolateSideAlongPolyline(
                t,
                l,
                startSD->n_plus,
                startSD->w_plus,
                blendEnd,
                endSD->n_minus,   // swap side at end (forward tangent is opposite to end outgoing)
                endSD->w_minus);
            const SideInterpolationResult minus = interpolateSideAlongPolyline(
                t,
                l,
                startSD->n_minus,
                startSD->w_minus,
                blendEnd,
                endSD->n_plus,    // swap side at end
                endSD->w_plus);
            plusN = plus.normals;
            minusN = minus.normals;
            plusW = plus.widths;
            minusW = minus.widths;
        } else {
            // Reverse interpolation from end intersection to start anchor, then map back.
            std::vector<Eigen::Vector3d> tRev(k, Eigen::Vector3d::Zero());
            std::vector<double> lRev(k, 0.0);
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t i = k - 1 - j;
                tRev[j] = -t[i];
                lRev[j] = l[i];
            }

            const bool blendEnd = false;  // one-intersection/anchor case
            const SideInterpolationResult plusRev = interpolateSideAlongPolyline(
                tRev, lRev, endSD->n_plus, endSD->w_plus, blendEnd, Eigen::Vector3d::Zero(), endSD->w_plus);
            const SideInterpolationResult minusRev = interpolateSideAlongPolyline(
                tRev, lRev, endSD->n_minus, endSD->w_minus, blendEnd, Eigen::Vector3d::Zero(), endSD->w_minus);

            for (std::size_t i = 0; i < k; ++i) {
                const std::size_t j = k - 1 - i;
                // Reverse orientation flips left/right assignment.
                plusN[i] = minusRev.normals[j];
                minusN[i] = plusRev.normals[j];
                plusW[i] = minusRev.widths[j];
                minusW[i] = plusRev.widths[j];
            }
        }

        for (std::size_t i = 0; i < k; ++i) {
            const int segIdx = segIdxs[i];
            if (segIdx < 0 || static_cast<std::size_t>(segIdx) >= dSegments.size()) {
                continue;
            }

            dsegment& seg = dSegments[static_cast<std::size_t>(segIdx)];
            const Eigen::Vector3d ti = seg.direction().normalized();
            const Eigen::Vector3d tPlus = ti;
            const Eigen::Vector3d tMinus = -ti;

            SegmentSideData plusSide;
            plusSide.n = plusN[i];
            plusSide.w = std::max(kWidthEps, plusW[i]);
            plusSide.h = std::sqrt(std::max(kWidthEps, seg.segmentLength() * plusSide.w));
            plusSide.b = plusSide.n.cross(tPlus).normalized();
            plusSide.B.col(0) = tPlus;
            plusSide.B.col(1) = plusSide.b;
            plusSide.B.col(2) = plusSide.n;
            plusSide.S = Eigen::Vector3d(seg.segmentLength(), plusSide.w, plusSide.h);
            plusSide.BS = plusSide.B * plusSide.S.asDiagonal();
            plusSide.valid = true;

            SegmentSideData minusSide;
            minusSide.n = minusN[i];
            minusSide.w = std::max(kWidthEps, minusW[i]);
            minusSide.h = std::sqrt(std::max(kWidthEps, seg.segmentLength() * minusSide.w));
            minusSide.b = minusSide.n.cross(tMinus).normalized();
            minusSide.B.col(0) = tMinus;
            minusSide.B.col(1) = minusSide.b;
            minusSide.B.col(2) = minusSide.n;
            minusSide.S = Eigen::Vector3d(seg.segmentLength(), minusSide.w, minusSide.h);
            minusSide.BS = minusSide.B * minusSide.S.asDiagonal();
            minusSide.valid = true;

            seg.setPlusSide(plusSide);
            seg.setMinusSide(minusSide);
        }
    }
}

void dcurvenet::initializeFromCurvenet(const Curvenet::curvenet& CN) {
    dVerts.clear();
    dSegments.clear();
    dSplines.clear();

    const auto& controls = CN.controls();
    const auto& splines = CN.getSplines();

    std::vector<int> controlToDvert(controls.size(), -1);
    dVerts.reserve(controls.size());

    for (std::size_t cIdx = 0; cIdx < controls.size(); ++cIdx) {
        const int dIdx = static_cast<int>(dVerts.size());
        controlToDvert[cIdx] = dIdx;
        dVerts.emplace_back(controls[cIdx].getPosition(), static_cast<int>(cIdx), -1, 0.0);
    }

    // Pass 1: copy each spline as a dspline placeholder, preserving index alignment.
    dSplines.reserve(splines.size());
    for (std::size_t sIdx = 0; sIdx < splines.size(); ++sIdx) {
        const Curvenet::spline& s = splines[sIdx];
        const int startControl = s.startControl();
        const int endControl = s.endControl();
        if (startControl < 0 || endControl < 0 ||
            static_cast<std::size_t>(startControl) >= controls.size() ||
            static_cast<std::size_t>(endControl) >= controls.size()) {
            throw std::runtime_error("Spline endpoint control index out of range");
        }

        const int startDvert = controlToDvert[static_cast<std::size_t>(startControl)];
        const int endDvert = controlToDvert[static_cast<std::size_t>(endControl)];
        dSplines.emplace_back(static_cast<int>(sIdx), startDvert, endDvert);
    }

    // Copy control->spline adjacency order from curvenet controls (CCW order as stored there).
    for (std::size_t cIdx = 0; cIdx < controls.size(); ++cIdx) {
        const int dIdx = controlToDvert[cIdx];
        if (dIdx < 0) {
            continue;
        }
        dVerts[static_cast<std::size_t>(dIdx)].setDSplineIdxs(controls[cIdx].getSplineIdxs());
    }

    // Pass 2: discretize each dspline and populate its segments/range pointers.
    for (std::size_t sIdx = 0; sIdx < splines.size(); ++sIdx) {
        const Curvenet::spline& s = splines[sIdx];
        dspline& ds = dSplines[sIdx];

        const int interiorCount = computeInteriorSampleCount(s, mean_edge_length, samples_per_mean_edge);
        const int outputCount = interiorCount + 2;
        const int refineCount = std::max(uniform_refine_samples, outputCount);

        const std::vector<Eigen::Vector3d> refined = s.sampleUniformParameterization(refineCount);
        const std::vector<Eigen::Vector3d> sampled = resamplePolylineEvenArcLength(refined, outputCount);

        std::vector<int> splineDverts;
        splineDverts.reserve(static_cast<std::size_t>(outputCount));
        splineDverts.push_back(ds.startDvert());

        for (int i = 1; i < outputCount - 1; ++i) {
            const int dIdx = static_cast<int>(dVerts.size());
            const double t = static_cast<double>(i) / static_cast<double>(outputCount - 1);
            dVerts.emplace_back(sampled[static_cast<std::size_t>(i)], -1, static_cast<int>(sIdx), t);
            splineDverts.push_back(dIdx);
        }

        splineDverts.push_back(ds.endDvert());

        std::vector<int> segmentIdxs;
        segmentIdxs.reserve(splineDverts.size() - 1);

        for (std::size_t i = 1; i < splineDverts.size(); ++i) {
            const int sDvert = splineDverts[i - 1];
            const int eDvert = splineDverts[i];

            const int segIdx = static_cast<int>(dSegments.size());
            dSegments.emplace_back(sDvert,
                                   eDvert,
                                   static_cast<int>(sIdx),
                                   dVerts[static_cast<std::size_t>(sDvert)].position(),
                                   dVerts[static_cast<std::size_t>(eDvert)].position());

            dVerts[static_cast<std::size_t>(sDvert)].addSegmentIdx(segIdx);
            dVerts[static_cast<std::size_t>(eDvert)].addSegmentIdx(segIdx);
            segmentIdxs.push_back(segIdx);
        }

        ds.setInteriorSamples(interiorCount);
        ds.setSegments(std::move(segmentIdxs));
    }

    computeControlRibbonData(controls, controlToDvert);
    interpolateRibbonAlongDSplines();
}

}  // namespace DCurvenet
