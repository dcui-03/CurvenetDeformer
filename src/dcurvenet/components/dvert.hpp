// dvert.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

namespace DCurvenet {

// Control Points
class dvert {
    public:
        // Constructor which takes the 
        dvert(Eigen::Vector3d position, int control);

        // TODO: Needs getters that check the labels and can return an error if we try to get something with wrong attribute

        // Set the normal (during parallel transport)
        // sign of true means positive
        int setNormal(Eigen::Vector3d normal, bool sign);
        // Number of outgoing segments from this control
        // Returns -1 if not a control
        int num_outgoing();
    protected:
        // No class inheritance
    private:
        // Attributes
        Eigen::Vector3d pos;
        int control;            // -1 if not a control, and the index of the parent control point if is
        std::vector<int> segment_idxs;   // pointers to an ordering of outgoing segments (CCW) or (previous then next along spline)
        // matching CCW corner normals for outgoing segments (if is from control)
        // For non-control vertices this should have positive first, then negative normal
        std::vector<Eigen::Vector3d> adjacent_normals;
        double t;       // t-value along parent spline.  NOTE: for control points, this must be 0.0 or 1.0
        
        // PDC ATTRIBUTES
        int DC_origin = 0;      // 0 for dverts inherited from DC, 1 for dverts produced by splits
        int projection_element; // 0 if projected dvert lands on face, 1 if lands on edge, 2 if lands on vertex
        int element_index;      // The index of the face, edge or vertex in the mesh
};

}   // namespace Curvenet