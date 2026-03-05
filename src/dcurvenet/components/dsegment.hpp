// dvert.hpp
#pragma once

#include "dvert.hpp"
#include "dspline.hpp"
#include <Eigen/Core>
#include <vector>

namespace DCurvenet {

// Control Points
class dsegment {
    public:
        // Constructor which takes the start and ends of the segment, as well as the parent dspline index
        dsegment(int start_idx, int end_idx, int dspline_idx);

        // TODO: Needs getters that check the labels and can return an error if we try to get something with wrong attribute

    protected:
        // No class inheritance
    private:
        // Compute Local Scaled Frames
        void computeLocalScalesandFrames();

        // Attributes
        Eigen::Vector3d direc;  // oriented direction of segment
        double length;  // length of segment
        // Local frames and their scales
        // Note, store scales as vector and turn into diagonal matrix when needed
        Eigen::Matrix3d p_frame;
        Eigen::Matrix3d m_frame;
        Eigen::Vector3d p_scale;
        Eigen::Vector3d m_scale;

        int start_dvert;
        int end_dvert;
        int parent_dspline;  // pointer to parent discrete spline

        // PDC ATTRIBUTES
        // Note: In PDC, we need to alter our start_dvert and end_dvert to be PDC's indices
        // Attributes we only need for PDC
        int DCsegment = -1; // Index of parent dsegment in DC
        
};

}   // namespace Curvenet