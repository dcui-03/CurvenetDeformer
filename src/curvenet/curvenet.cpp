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

curvenet::curvenet(const std::vector<Eigen::Vector3d>& controlP,
                   const std::vector<Eigen::Vector3d>& controlNormals,
                   const std::vector<Eigen::Vector3d>& tangentP,
                   const std::vector<std::array<int, 4>>& curveC) {
    initializeControls(controlP, controlNormals);

    for (const auto& c : curveC) {
        if (c[0] < 0 || c[3] < 0 ||
            static_cast<std::size_t>(c[0]) >= controlPoints.size() ||
            static_cast<std::size_t>(c[3]) >= controlPoints.size() ||
            c[1] < 0 || c[2] < 0 ||
            static_cast<std::size_t>(c[1]) >= tangentP.size() ||
            static_cast<std::size_t>(c[2]) >= tangentP.size()) {
            throw std::out_of_range("curve constructor indices out of range");
        }
        addSpline(std::array<int, 2>{c[0], c[3]},
                  std::array<Eigen::Vector3d, 2>{tangentP[static_cast<std::size_t>(c[1])],
                                                 tangentP[static_cast<std::size_t>(c[2])]});
    }
}

control& curvenet::controlAt(int idx) {
    if (idx < 0 || static_cast<std::size_t>(idx) >= controlPoints.size()) {
        throw std::out_of_range("control index out of range");
    }
    return controlPoints[static_cast<std::size_t>(idx)];
}

const control& curvenet::controlAt(int idx) const {
    if (idx < 0 || static_cast<std::size_t>(idx) >= controlPoints.size()) {
        throw std::out_of_range("control index out of range");
    }
    return controlPoints[static_cast<std::size_t>(idx)];
}

tangent& curvenet::tangentAt(int idx) {
    if (idx < 0 || static_cast<std::size_t>(idx) >= tangentPoints.size()) {
        throw std::out_of_range("tangent index out of range");
    }
    return tangentPoints[static_cast<std::size_t>(idx)];
}

const tangent& curvenet::tangentAt(int idx) const {
    if (idx < 0 || static_cast<std::size_t>(idx) >= tangentPoints.size()) {
        throw std::out_of_range("tangent index out of range");
    }
    return tangentPoints[static_cast<std::size_t>(idx)];
}

spline& curvenet::splineAt(int idx) {
    if (idx < 0 || static_cast<std::size_t>(idx) >= splines.size()) {
        throw std::out_of_range("spline index out of range");
    }
    return splines[static_cast<std::size_t>(idx)];
}

const spline& curvenet::splineAt(int idx) const {
    if (idx < 0 || static_cast<std::size_t>(idx) >= splines.size()) {
        throw std::out_of_range("spline index out of range");
    }
    return splines[static_cast<std::size_t>(idx)];
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
        sourceToControl[i] = addControl(controlP[i], n);
    }
}

void curvenet::initializeControls(const std::vector<Eigen::Vector3d>& controlP,
                                  const std::vector<Eigen::Vector3d>& controlNormals) {
    controlPoints.clear();
    tangentPoints.clear();
    splines.clear();
    sourceToControl.clear();

    if (controlP.size() != controlNormals.size()) {
        throw std::invalid_argument("controlP and controlNormals must have matching sizes");
    }
    for (std::size_t i = 0; i < controlP.size(); ++i) {
        addControl(controlP[i], controlNormals[i]);
    }
}

void curvenet::intializeSplines(const std::vector<Eigen::Vector3d>& controlP,
                                const std::vector<std::vector<int>>& curveC) {
    for (const auto& c : curveC) {
        const int c1 = controlIndexFromSource(c[0]);
        const int c2 = controlIndexFromSource(c[3]);
        if (c1 == c2) {
            continue;
        }

        addSpline(std::array<int, 2>{c1, c2},
                  std::array<Eigen::Vector3d, 2>{controlP[static_cast<std::size_t>(c[1])],
                                                 controlP[static_cast<std::size_t>(c[2])]});
    }
}

int curvenet::addControl(const Eigen::Vector3d& controlP, const Eigen::Vector3d& controlNormal) {
    controlPoints.emplace_back(controlP, controlNormal);
    return static_cast<int>(controlPoints.size()) - 1;
}

int curvenet::appendSpline(int control0, int tangent0, int tangent1, int control1) {
    const int splineIdx = static_cast<int>(splines.size());
    tangentPoints[static_cast<std::size_t>(tangent0)].setSplineIdx(splineIdx);
    tangentPoints[static_cast<std::size_t>(tangent1)].setSplineIdx(splineIdx);

    splines.emplace_back(std::vector<int>{control0, tangent0, tangent1, control1});
    splines.back().computeSpline(controlPoints, tangentPoints);

    controlPoints[static_cast<std::size_t>(control0)].addSplineIdx(splineIdx);
    controlPoints[static_cast<std::size_t>(control1)].addSplineIdx(splineIdx);
    reorderControlSplineOrderingCCW(control0);
    if (control1 != control0) {
        reorderControlSplineOrderingCCW(control1);
    }
    return splineIdx;
}

int curvenet::addSpline(const std::array<Eigen::Vector3d, 2>& controlP,
                        const std::array<Eigen::Vector3d, 2>& controlNormals,
                        const std::array<Eigen::Vector3d, 2>& tangentP) {
    const int control0 = addControl(controlP[0], controlNormals[0]);
    const int control1 = addControl(controlP[1], controlNormals[1]);
    return addSpline(std::array<int, 2>{control0, control1}, tangentP);
}

int curvenet::addSpline(const std::array<int, 2>& controlIdx,
                        const std::array<Eigen::Vector3d, 2>& tangentP) {
    const int control0 = controlIdx[0];
    const int control1 = controlIdx[1];
    if (control0 < 0 || control1 < 0 ||
        static_cast<std::size_t>(control0) >= controlPoints.size() ||
        static_cast<std::size_t>(control1) >= controlPoints.size()) {
        throw std::out_of_range("control index out of range");
    }

    const int tangent0 = static_cast<int>(tangentPoints.size());
    tangentPoints.emplace_back(tangentP[0], control0);

    const int tangent1 = static_cast<int>(tangentPoints.size());
    tangentPoints.emplace_back(tangentP[1], control1);

    return appendSpline(control0, tangent0, tangent1, control1);
}

int curvenet::addSpline(int control0,
                        const Eigen::Vector3d& control1,
                        const Eigen::Vector3d& normal1,
                        const std::array<Eigen::Vector3d, 2>& tangentP) {
    const int newControl = addControl(control1, normal1);
    return addSpline(std::array<int, 2>{control0, newControl}, tangentP);
}

int curvenet::addSpline(const Eigen::Vector3d& control0,
                        int control1,
                        const Eigen::Vector3d& normal0,
                        const std::array<Eigen::Vector3d, 2>& tangentP) {
    const int newControl = addControl(control0, normal0);
    return addSpline(std::array<int, 2>{newControl, control1}, tangentP);
}

void curvenet::setControlPosition(int controlIdx, const Eigen::Vector3d& p) {
    controlAt(controlIdx).setPosition(p);
    recomputeAllSplines();
}

void curvenet::setTangentPosition(int tangentIdx, const Eigen::Vector3d& p) {
    tangent& t = tangentAt(tangentIdx);
    t.setPosition(p);
    const int splineIdx = t.getSplineIdx();
    if (splineIdx >= 0) {
        recomputeSpline(splineIdx);
    } else {
        recomputeAllSplines();
    }
}

void curvenet::recomputeSpline(int splineIdx) {
    splineAt(splineIdx).computeSpline(controlPoints, tangentPoints);
    const spline& s = splineAt(splineIdx);
    reorderControlSplineOrderingCCW(s.startControl());
    reorderControlSplineOrderingCCW(s.endControl());
}

void curvenet::recomputeAllSplines() {
    for (auto& s : splines) {
        s.computeSpline(controlPoints, tangentPoints);
    }
    reorderControlSplineOrderingCCW();
}

void curvenet::reorderControlSplineOrderingCCW() {
    for (std::size_t cIdx = 0; cIdx < controlPoints.size(); ++cIdx) {
        reorderControlSplineOrderingCCW(static_cast<int>(cIdx));
    }
}

void curvenet::reorderControlSplineOrderingCCW(int controlIdx) {
    constexpr double kEps = 1e-12;

    if (controlIdx < 0 || static_cast<std::size_t>(controlIdx) >= controlPoints.size()) {
        return;
    }

    const std::vector<int> current = controlPoints[static_cast<std::size_t>(controlIdx)].getSplineIdxs();
    if (current.size() <= 1) {
        return;
    }

    const Eigen::Vector3d n = controlPoints[static_cast<std::size_t>(controlIdx)].getNormal().normalized();

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

        Eigen::Vector3d t = splines[static_cast<std::size_t>(sIdx)].outgoingTangentAtControl(controlIdx);
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
        return;
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
    controlPoints[static_cast<std::size_t>(controlIdx)].setSplineIdxs(sorted);
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

void curvenet::convertControlsToPC(std::vector<Eigen::Vector3d>& psControls) const {
    psControls.clear();
    psControls.reserve(controlPoints.size());
    for (const auto& c : controlPoints) {
        psControls.push_back(c.getPosition());
    }
}

void curvenet::convertTangentsToPC(std::vector<Eigen::Vector3d>& psTangents) const {
    psTangents.clear();
    psTangents.reserve(tangentPoints.size());
    for (const auto& t : tangentPoints) {
        psTangents.push_back(t.getPosition());
    }
}

void curvenet::convertCurvnetToCN(std::vector<Eigen::Vector3d>& psSamples,
                                  std::vector<std::array<int, 2>>& psE,
                                  int samplesPerSpline) const {
    psE.clear();
    convertControlsToPC(psSamples);

    for (const auto& s : splines) {
        const std::vector<Eigen::Vector3d> samples = s.sampleParameterization(std::max(2, samplesPerSpline));
        int prevIdx = s.startControl();
        for (std::size_t i = 1; i + 1 < samples.size(); ++i) {
            const int vIdx = static_cast<int>(psSamples.size());
            psSamples.push_back(samples[i]);
            psE.push_back({prevIdx, vIdx});
            prevIdx = vIdx;
        }
        psE.push_back({prevIdx, s.endControl()});
    }
}

void curvenet::convertControlsAndTangentsToCN(std::vector<Eigen::Vector3d>& psCAndT,
                                              std::vector<std::array<int, 2>>& psCAndTE) const {
    psCAndTE.clear();
    convertControlsToPC(psCAndT);
    for (const auto& s : splines) {
        const int startTanIdx = static_cast<int>(psCAndT.size());
        psCAndT.push_back(tangentPoints[static_cast<std::size_t>(s.startTangent())].getPosition());
        psCAndTE.push_back({s.startControl(), startTanIdx});

        const int endTanIdx = static_cast<int>(psCAndT.size());
        psCAndT.push_back(tangentPoints[static_cast<std::size_t>(s.endTangent())].getPosition());
        psCAndTE.push_back({s.endControl(), endTanIdx});
    }
}

void curvenet::modifyPositions(const std::vector<Eigen::Vector3d>& controlP) {
    if (controlP.size() != controlPoints.size()) {
        return;
    }
    for (std::size_t i = 0; i < controlPoints.size(); ++i) {
        controlPoints[i].setPosition(controlP[i]);
    }
    recomputeAllSplines();
}

void curvenet::modifyPositions(const std::vector<Eigen::Vector3d>& newControls,
                               const std::vector<Eigen::Vector3d>& newTangents) {
    if (newControls.size() == controlPoints.size()) {
        for (std::size_t i = 0; i < controlPoints.size(); ++i) {
            controlPoints[i].setPosition(newControls[i]);
        }
    }
    if (newTangents.size() == tangentPoints.size()) {
        for (std::size_t i = 0; i < tangentPoints.size(); ++i) {
            tangentPoints[i].setPosition(newTangents[i]);
        }
    }
    recomputeAllSplines();
}

}  // namespace Curvenet
