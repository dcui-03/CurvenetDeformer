// curvenet.hpp
#pragma once

#include "components/control.hpp"
#include "components/spline.hpp"

#include <Eigen/Core>

#include <vector>


namespace Curvenet {

class curvenet {
public:
    // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
    // TODO change structure of controls?
    curvenet(const std::vector<Eigen::Vector3d>& controlP,
             const std::vector<Eigen::Vector3d>& surfaceN,
             const std::vector<std::vector<int>>& curveC);

    // TODO: Needs getters so that the discrete CN can ask for internals
    const std::vector<control>& controls() const { return controlPoints; }
    const std::vector<tangent>& tangents() const { return tangentPoints; }
    const std::vector<spline>& getSplines() const { return splines; }

    // Legacy stub signatures preserved from original header:
    // modifyPositions();
    // initializePoints();
    // intializeSplines();

    // Modify the positions of the curvenet for deforming
    void modifyPositions(const std::vector<Eigen::Vector3d>& controlP);

protected:
    // No class inheritance
private:
    // Reorganize during intialization
    void initializePoints(const std::vector<Eigen::Vector3d>& controlP,
                          const std::vector<Eigen::Vector3d>& surfaceN,
                          const std::vector<std::vector<int>>& curveC);
    void intializeSplines(const std::vector<Eigen::Vector3d>& controlP,
                          const std::vector<std::vector<int>>& curveC);
    int controlIndexFromSource(int sourceIdx) const;

    // Store control points as a list
    std::vector<control> controlPoints;
    std::vector<tangent> tangentPoints;
    std::vector<spline> splines;

    // Mapping from source point indices to internally stored controls (endpoints only)
    std::vector<int> sourceToControl;
};

}  // namespace Curvenet
