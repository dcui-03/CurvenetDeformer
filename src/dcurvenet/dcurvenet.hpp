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
        // NOTE: Vertices are copied directly from the Control of the curve network (CN),
        //       meaning they are in the same order and have the same corresponding indices.
        //       Curves are similar.
        // Takes the original curvenet and discretizes it
        // Alpha is the user-inputted sampling parameter
        dcurvenet(Curvenet::curvenet* CN, int alpha = 5);
        // Initialize with empty constructor
        dcurvenet();

        // Compute the deformation gradient on an edge given a new scaled frame
        // NOTE: Use formula from paper
        void computeHEDefGrad(int he, Eigen::Vector3d newT,
                                      Eigen::Vector3d newB,
                                      Eigen::Vector3d newN,
                                      Eigen::Vector3d newScale);

    protected:
        // No inherited classes
    private:

        // For controls, computes their corner normals and widths. For non-intersections, this method does nothing (return -1)
        int vertCornerNormalsWidths(int v);
        // Corner normals on all vertices
        int allCornerNormalsAndWidths();
        // Transport corner normals and widths from the two end corners of a curve
        int transportNWOnCurve(int c);
        // Transport normals and widths for all curves
        int transportNormalsAndWidths();

        // Compute local frame on a curve
        int computeScaledFrameOnCurve(int c);   // Notice that we need to do this
        // Compute scaled frames on all curves
        int computeScaledFrames();
            
        // Attributes as lists
        std::vector<Vert> V;
        std::vector<Edge> E;
        std::vector<HalfEdge> HE;
        std::vector<Curve> C;

        // Pointer to parent curvenet object
        Curvenet::curvenet* CN;
};

}   // namespace DCurvenet