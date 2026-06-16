// dcurvenet.hpp
#pragma once

#include "dcurvenet_types.hpp"
#include "curvenet/curvenet.hpp"
#include <Eigen/Core>
#include <vector>

namespace DCurvenet {

// Discrete curvenet (i.e., polylines)
class dcurvenet {
    public:
        // Takes the original curvenet and discretizes it
        // Alpha is the user-inputted sampling parameter
        dcurvenet(const Curvenet::curvenet& CN, int alpha = 5);
        // Initialize with empty constructor
        dcurvenet();

        // Compute the deformation gradient on an edge given a new scaled frame
        // NOTE: Use formula from paper
        void computeHEDefGrad(int he, Eigen::Matrix3d newFrame, Eigen::Vector3d newScale);

    protected:
        // No inherited classes
    private:

        // For controls, computes their corner normals. For non-controls, this method does nothing (return -1)
        int vertCornerNormals(int v);
        // Corner normals on all vertices
        int allCornerNormals();
        // Transport corner normals from the two ends of a curve
        int transportNormalOnCurve(int c);
        // Transport normals for all curves
        int transportNormals();

        // Compute local frame on a curve
        int computeScaledFrameOnCurve(int c);   // Notice that we need to do this
        // Compute scaled frames on all curves
        int computeScaledFrames();
        
        std::vector<Vert> V;
        std::vector<Edge> E;
        std::vector<HalfEdge> HE;
        std::vector<Curve> C;
};

}   // namespace DCurvenet