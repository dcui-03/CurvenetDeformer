// curvenet.hpp
#pragma once

#include "components/control.hpp"
#include "components/spline.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace Curvenet {

class curvenet {
    public:
        // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
        // TODO change structure of controls?
        curvenet(std::vector<Eigen::Vector3d> controlP, std::vector<Eigen::Vector3d> surfaceN, std::vector<std::vector<int>> curveC);

        // TODO: Needs getters so that the discrete CN can ask for internals

        // Modify the positions of the curvenet for deforming
        modifyPositions();
    protected:
        // No class inheritance
    private:
        // Reorganize during intialization
        initializePoints();
        intializeSplines();


        // Store control points as a list
        std::vector<control> controlPoints;
        std::vector<tangent> tangentPoints;
        std::vector<spline> splines;
        
        
};

}   // namespace CutMesh