// dspline.hpp
#pragma once

#include "dsegment.hpp"
#include "dvert.hpp"
#include <Eigen/Core>
#include <vector>


namespace DCurvenet {

class dspline {
    public:
        // constructor, which takes a list of segments
        dspline(std::vector<int> segments, int spline_idx);
    protected:
        // No class inheritance
    private:
        // Propagate normals between two endpoints (may need several variants to catch all cases)
        void parallelTransportNormals(Eigen::Vector3d startN, Eigen::Vector3d endN);

        // Store pointers to all segments in the spline
        std::vector<int> segments;
        // Also pointers to the control dverts (so it's easy to query)
        int start_dvert;
        int end_dvert;
        // Number of dvert samples
        int n;
        // Pointer to parent spline
        int parent_spline;
        
        
};

}   // namespace CutMesh