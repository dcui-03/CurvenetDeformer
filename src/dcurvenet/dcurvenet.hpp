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
        // NOTE: Uses formula from paper
        Eigen::Matrix3d computeHEDefGrad(int he,
                                        const Eigen::Vector3d& newT, 
                                        const Eigen::Vector3d& newB, 
                                        const Eigen::Vector3d& newN, 
                                        const double& newL, const double& newW, const double& newH);

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
        
        // Accumulates rotation matrices and lengths by tracing from a starting halfedge to an end vertex
        double accumulateRotations(int start_he, int end_v, std::vector<Eigen::Matrix3d>& rots, std::vector<double> lens);
        // Compute torsioin
        double computeTorsion(Eigen::Vector3d n_1, Eigen::Vector3d n_k, Eigen::Matrix3d Om_k, Eigen::Vector3d t_k);
            
        // Attributes as lists
        std::vector<Vert> V;
        std::vector<Edge> E;
        std::vector<HalfEdge> HE;
        std::vector<Curve> C;

        // Pointer to parent curvenet object
        Curvenet::curvenet* CN;
};

}   // namespace DCurvenet