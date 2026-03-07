#include "curvenet.hpp"

#include <limits>
#include <stdexcept>


namespace Curvenet {

namespace {
Eigen::Vector3d safeNormalize(const Eigen::Vector3d& v, const Eigen::Vector3d& fallback) {
    const double n = v.norm();
    if (n <= std::numeric_limits<double>::epsilon()) {
        return fallback;
    }
    return v / n;
}
}  // namespace

curvenet::curvenet(const std::vector<Eigen::Vector3d>& controlP,
                   const std::vector<Eigen::Vector3d>& surfaceN,
                   const std::vector<std::vector<int>>& curveC) {
    initializePoints(controlP, surfaceN, curveC);
    intializeSplines(controlP, curveC);
}

void curvenet::initializePoints(const std::vector<Eigen::Vector3d>& controlP,
                                const std::vector<Eigen::Vector3d>& surfaceN,
                                const std::vector<std::vector<int>>& curveC) {
    controlPoints.clear();
    tangentPoints.clear();
    splines.clear();
    sourceToControl.assign(controlP.size(), -1);

    std::vector<char> isEndpoint(controlP.size(), 0);
    for (const auto& c : curveC) {
        if (c.size() != 4) {
            throw std::invalid_argument("curveC rows must have 4 indices [start,t1,t2,end]");
        }
        for (int idx : c) {
            if (idx < 0 || static_cast<std::size_t>(idx) >= controlP.size()) {
                throw std::out_of_range("curveC index out of range");
            }
        }
        isEndpoint[static_cast<std::size_t>(c[0])] = 1;
        isEndpoint[static_cast<std::size_t>(c[3])] = 1;
    }

    for (std::size_t i = 0; i < controlP.size(); ++i) {
        if (!isEndpoint[i]) {
            continue;
        }
        Eigen::Vector3d n = Eigen::Vector3d::UnitZ();
        if (i < surfaceN.size()) {
            n = safeNormalize(surfaceN[i], n);
        } else {
            n = safeNormalize(controlP[i], n);
        }
        sourceToControl[i] = static_cast<int>(controlPoints.size());
        controlPoints.emplace_back(controlP[i], n);
    }
}

void curvenet::intializeSplines(const std::vector<Eigen::Vector3d>& controlP,
                                const std::vector<std::vector<int>>& curveC) {
    for (const auto& c : curveC) {
        const int c1 = controlIndexFromSource(c[0]);
        const int c2 = controlIndexFromSource(c[3]);
        if (c1 == c2) {
            continue;  // skip degenerate loops
        }

        const int t1 = static_cast<int>(tangentPoints.size());
        tangentPoints.emplace_back(controlP[static_cast<std::size_t>(c[1])], c1);

        const int t2 = static_cast<int>(tangentPoints.size());
        tangentPoints.emplace_back(controlP[static_cast<std::size_t>(c[2])], c2);

        const int sIdx = static_cast<int>(splines.size());
        splines.emplace_back(std::vector<int>{c1, t1, t2, c2});
        splines.back().computeSpline(controlPoints, tangentPoints);

        tangentPoints[static_cast<std::size_t>(t1)].setSplineIdx(sIdx);
        tangentPoints[static_cast<std::size_t>(t2)].setSplineIdx(sIdx);

        controlPoints[static_cast<std::size_t>(c1)].addSplineIdx(sIdx);
        controlPoints[static_cast<std::size_t>(c2)].addSplineIdx(sIdx);
    }
}

int curvenet::controlIndexFromSource(int sourceIdx) const {
    if (sourceIdx < 0 || static_cast<std::size_t>(sourceIdx) >= sourceToControl.size()) {
        throw std::out_of_range("source control index out of range");
    }
    const int mapped = sourceToControl[static_cast<std::size_t>(sourceIdx)];
    if (mapped < 0) {
        throw std::runtime_error("source endpoint does not map to a control");
    }
    return mapped;
}

void curvenet::modifyPositions(const std::vector<Eigen::Vector3d>& controlP) {
    // Not needed for curvenet initialization/render path yet.
    // Keeping stub to preserve intended API shape.
    if (controlP.size() != controlPoints.size()) {
        return;
    }
    for (std::size_t i = 0; i < controlPoints.size(); ++i) {
        controlPoints[i].setPosition(controlP[i]);
    }
    for (auto& s : splines) {
        s.computeSpline(controlPoints, tangentPoints);
    }
}

}  // namespace Curvenet
