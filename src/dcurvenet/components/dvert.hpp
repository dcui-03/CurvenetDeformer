// dvert.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <utility>

namespace DCurvenet {

// Control Points
class dvert {
    public:
        // Constructor which takes the position
        dvert(); 
        dvert(Eigen::Vector3d position);
        dvert(Eigen::Vector3d position, int origin, int element_type, int element_index);

        // TODO: Needs getters that check the labels and can return an error if we try to get something with wrong attribute
        // Set position after initialization
        int setPosition(Eigen::Vector3d position);
        Eigen::Vector3d dvert::getPosition();
        // Set the normal (during parallel transport)
        // sign of true means positive
        int setNormal(Eigen::Vector3d normal, bool sign);
        // Number of outgoing segments from this control
        // Returns -1 if not a control
        int num_outgoing();
        // Returns the projection element and element idx
        std::pair<int, int> getProjection();
        // Set the projection element
        void setProjection(int proj_element, int el_idx);
    protected:
        // Make a hard copy (ex. for PDC)
        dvert makeCopy();

        // Attributes
        Eigen::Vector3d pos;
        bool control = false;            // bool that is true if is control, and false otherwise
        std::vector<int> dspline_idxs;   // pointer to parent splines (CCW) or just the 1.
        std::vector<int> dsegment_idxs;   // pointers to an ordering of outgoing segments (CCW) or (previous then next along spline)
        // matching CCW corner normals for outgoing segments (if is from control)
        // For non-control vertices this should have positive first, then negative normal
        std::vector<Eigen::Vector3d> adjacent_normals;
        double t;       // t-value along parent spline.  NOTE: for control points, this must be 0.0 or 1.0
        
        // PDC ATTRIBUTES
        int DC_origin = 0;      // 0 for dverts inherited from DC, 1 for dverts produced by splits
        int projection_element; // 0 if projected dvert lands on face, 1 if lands on edge, 2 if lands on vertex
        int element_idx;      // The index of the face, edge or vertex in the mesh
    private:
};

/*
class dcontrol : public dvert {
    public:
        // Takes position and parent control index
        dcontrol(Eigen::Vector3d position, int control);
    protected:
    private:
}
*/
}   // namespace Curvenet