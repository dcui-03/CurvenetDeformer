// dcurvenet.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "components/dvert.hpp"
#include "components/dsegment.hpp"
#include "components/dspline.hpp"
#include <Eigen/Core>


namespace DCurvenet {

class dcurvenet {
    public:
        // Takes the original curvenet and discretizes it
        dcurvenet(Curvenet::curvenet& CN);
        // Initialize with empty constructor
        dcurvenet();

        // TODO: Needs getters so that others can ask for internals
        dvert getDVert(int idx);
        Eigen::Vector3d getDVertPos(int idx);
        dsegment getDSegment(int idx);
        dspline getDSpline(int idx);

    protected:
        // Reorganize during intialization
        void intializeDSplines(Curvenet::curvenet& CN);
        void initializeDSegments(Curvenet::curvenet& CN);
        void initializeDVerts(Curvenet::curvenet& CN);

        // Helper to set this as a copy of 

        // Helper to set this as a copy of another dcurvenet
        void copyFromDCurvenet(dcurvenet& source);

        // split a dsegment (by index) at a specified dvert
        void splitDSegment(dvert splitPoint, int splitSegment);

        // add a dvert to the list and return its index
        int addDVert(Eigen::Vector3d addPoint);


        // Store control points as a list
        // Store verts as first copying the curvenet. Then add middle points after that spline by spline
        // This preserves indexing for controls
        std::vector<dvert> dVerts;
        std::vector<dsegment> dSegments;
        std::vector<dspline> dSplines;
    private:
        // For controls, computes their corner normals. For non-controls, this method does nothing (return -1)
        int computeCornerNormals();
        
        
};

}   // namespace DCurvenet