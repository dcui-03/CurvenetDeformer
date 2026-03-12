// curvenet.hpp
#pragma once

#include "components/control.hpp"
#include "components/spline.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <map>
#include <array>

namespace Curvenet {

class curvenet {
    public:
        // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
        // NOTE: Constructor assumes you already have no duplicates in your inputs
        curvenet(std::vector<Eigen::Vector3d> controlP, 
                 std::vector<Eigen::Vector3d> controlNormals, 
                 std::vector<Eigen::Vector3d> tangentP,
                 std::vector<std::array<int, 4>> curveC);
        // Empty constructor
        curvenet();
        
        // Add control or spline. (Tangents must be added with a spline)
        void addControl(Eigen::Vector3d controlP, Eigen::Vector3d controlNormal);

        // Overloads for various situations
        void addSpline(std::array<Eigen::Vector3d, 2> controlP,
                       std::array<Eigen::Vector3d, 2> controlNormals,
                       std::array<Eigen::Vector3d, 2> tangentP);
        
        // Both controls already exist
        void addSpline(std::array<int, 2> controlP, std::array<Eigen::Vector3d, 2> tangentP);
        // First control exists
        void addSpline(int control0, Eigen::Vector3d control1, Eigen::Vector3d normal1, std::array<Eigen::Vector3d, 2> tangentP);
        // Second control exists
        void addSpline(Eigen::Vector3d control0, int control1, Eigen::Vector3d normal0, std::array<Eigen::Vector3d, 2> tangentP);
        
        // Convert a curvenet object into something polyscope can read
        // Try to preserve control/tangent indexing as much as possible
        void convertControlsToPC(std::vector<Eigen::Vector3d>& psControls);
        void convertTangentsToPC(std::vector<Eigen::Vector3d>& psTangents);
        void convertCurvnetToCN(std::vector<Eigen::Vector3d>& psCAndT, std::vector<std::array<int, 2>>& psE);
        void convertControlsAndTangentsToCN(std::vector<Eigen::Vector3d>& psCAndT, std::vector<std::array<int, 2>>& psCAndTE);

        // Modify the positions of the curvenet for deforming
        void moveControl(); // Don't forget to modify tangents WITH controls

        void modifyPositions(std::vector<Eigen::Vector3d> newControls, std::vector<Eigen::Vector3d> newTangents);


        // Store control points as a list
        std::vector<control> controlPoints;
        std::vector<tangent> tangentPoints;
        std::vector<spline> splines;
    protected:
        // No class inheritance
    private:
        // Add controls 
        void initializeControls(std::vector<Eigen::Vector3d>& controlP,
                                std::vector<Eigen::Vector3d>& controlNormals);

        // Compute CCW spline ordering for each control or all
        void computeCCWOrderingAll();
        void computeCCWOrdering(int ctrl_idx);
        
        double tol = 1e-5;
        
};

}   // namespace CutMesh