#include "curvenet.hpp"
#include "utils/utils.hpp"
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace Curvenet {

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
        const Eigen::Vector3d n = (i < surfaceN.size() ? surfaceN[i] : controlP[i]).normalized();
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

    reorderControlSplineOrderingCCW();
}

void curvenet::reorderControlSplineOrderingCCW() {
    constexpr double kEps = 1e-12;

    for (std::size_t cIdx = 0; cIdx < controlPoints.size(); ++cIdx) {
        const std::vector<int> current = controlPoints[cIdx].getSplineIdxs();
        if (current.size() <= 1) {
            continue;
        }

        const Eigen::Vector3d n = controlPoints[cIdx].getNormal().normalized();

        struct Entry {
            int splineIdx = -1;
            Eigen::Vector3d tProj = Eigen::Vector3d::Zero();
            double angle = 0.0;
            bool degenerate = false;
            std::size_t stableIdx = 0;
        };

        std::vector<Entry> entries;
        entries.reserve(current.size());

        Eigen::Vector3d xRef = Eigen::Vector3d::Zero();
        bool haveXRef = false;

        for (std::size_t i = 0; i < current.size(); ++i) {
            const int sIdx = current[i];
            if (sIdx < 0 || static_cast<std::size_t>(sIdx) >= splines.size()) {
                continue;
            }

            Eigen::Vector3d t = splines[static_cast<std::size_t>(sIdx)].outgoingTangentAtControl(static_cast<int>(cIdx));
            Eigen::Vector3d tProj = t - n * t.dot(n);
            const double tn = tProj.norm();
            bool degenerate = tn <= kEps;
            if (!degenerate) {
                tProj /= tn;
                if (!haveXRef) {
                    xRef = tProj;
                    haveXRef = true;
                }
            } else {
                tProj = Eigen::Vector3d::Zero();
            }

            entries.push_back({sIdx, tProj, 0.0, degenerate, i});
        }

        if (entries.size() <= 1) {
            continue;
        }

        if (!haveXRef) {
            xRef = Utils::anyUnitTangent(n);
        }
        const Eigen::Vector3d yRef = n.cross(xRef).normalized();
        xRef = yRef.cross(n).normalized();

        for (auto& e : entries) {
            if (e.degenerate) {
                e.angle = std::numeric_limits<double>::infinity();
                continue;
            }
            e.angle = std::atan2(e.tProj.dot(yRef), e.tProj.dot(xRef));
        }

        std::stable_sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
            if (a.degenerate != b.degenerate) {
                return !a.degenerate;
            }
            if (a.angle == b.angle) {
                return a.stableIdx < b.stableIdx;
            }
            return a.angle < b.angle;
        });

        std::vector<int> sorted;
        sorted.reserve(entries.size());
        for (const Entry& e : entries) {
            sorted.push_back(e.splineIdx);
        }
        controlPoints[cIdx].setSplineIdxs(sorted);
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
