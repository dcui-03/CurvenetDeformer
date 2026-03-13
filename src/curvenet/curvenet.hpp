// curvenet.hpp
#pragma once

#include "components/control.hpp"
#include "components/spline.hpp"

#include <Eigen/Core>

#include <array>
#include <vector>

namespace Curvenet {

class curvenet {
public:
    curvenet() = default;
    curvenet(const std::vector<Eigen::Vector3d>& controlP,
             const std::vector<Eigen::Vector3d>& surfaceN,
             const std::vector<std::vector<int>>& curveC);
    curvenet(const std::vector<Eigen::Vector3d>& controlP,
             const std::vector<Eigen::Vector3d>& controlNormals,
             const std::vector<Eigen::Vector3d>& tangentP,
             const std::vector<std::array<int, 4>>& curveC);

    const std::vector<control>& controls() const { return controlPoints; }
    std::vector<control>& controlsMutable() { return controlPoints; }
    const std::vector<tangent>& tangents() const { return tangentPoints; }
    std::vector<tangent>& tangentsMutable() { return tangentPoints; }
    const std::vector<spline>& getSplines() const { return splines; }
    std::vector<spline>& splinesMutable() { return splines; }

    control& controlAt(int idx);
    const control& controlAt(int idx) const;
    tangent& tangentAt(int idx);
    const tangent& tangentAt(int idx) const;
    spline& splineAt(int idx);
    const spline& splineAt(int idx) const;

    int addControl(const Eigen::Vector3d& controlP, const Eigen::Vector3d& controlNormal);
    int addSpline(const std::array<Eigen::Vector3d, 2>& controlP,
                  const std::array<Eigen::Vector3d, 2>& controlNormals,
                  const std::array<Eigen::Vector3d, 2>& tangentP);
    int addSpline(const std::array<int, 2>& controlIdx,
                  const std::array<Eigen::Vector3d, 2>& tangentP);
    int addSpline(int control0,
                  const Eigen::Vector3d& control1,
                  const Eigen::Vector3d& normal1,
                  const std::array<Eigen::Vector3d, 2>& tangentP);
    int addSpline(const Eigen::Vector3d& control0,
                  int control1,
                  const Eigen::Vector3d& normal0,
                  const std::array<Eigen::Vector3d, 2>& tangentP);

    void setControlPosition(int controlIdx, const Eigen::Vector3d& p);
    void setTangentPosition(int tangentIdx, const Eigen::Vector3d& p);
    void recomputeSpline(int splineIdx);
    void recomputeAllSplines();

    void convertControlsToPC(std::vector<Eigen::Vector3d>& psControls) const;
    void convertTangentsToPC(std::vector<Eigen::Vector3d>& psTangents) const;
    void convertCurvnetToCN(std::vector<Eigen::Vector3d>& psSamples,
                            std::vector<std::array<int, 2>>& psE,
                            int samplesPerSpline = 64) const;
    void convertControlsAndTangentsToCN(std::vector<Eigen::Vector3d>& psCAndT,
                                        std::vector<std::array<int, 2>>& psCAndTE) const;

    void modifyPositions(const std::vector<Eigen::Vector3d>& controlP);
    void modifyPositions(const std::vector<Eigen::Vector3d>& newControls,
                         const std::vector<Eigen::Vector3d>& newTangents);

private:
    void initializePoints(const std::vector<Eigen::Vector3d>& controlP,
                          const std::vector<Eigen::Vector3d>& surfaceN,
                          const std::vector<std::vector<int>>& curveC);
    void initializeControls(const std::vector<Eigen::Vector3d>& controlP,
                            const std::vector<Eigen::Vector3d>& controlNormals);
    void intializeSplines(const std::vector<Eigen::Vector3d>& controlP,
                          const std::vector<std::vector<int>>& curveC);

    void reorderControlSplineOrderingCCW();
    void reorderControlSplineOrderingCCW(int controlIdx);
    int controlIndexFromSource(int sourceIdx) const;
    int appendSpline(int control0, int tangent0, int tangent1, int control1);

    std::vector<control> controlPoints;
    std::vector<tangent> tangentPoints;
    std::vector<spline> splines;
    std::vector<int> sourceToControl;
};

}  // namespace Curvenet
