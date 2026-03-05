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

        // TODO: Needs getters so that others can ask for internals

    protected:
        // No class inheritance
    private:
        // Reorganize during intialization
        intializeDSplines();
        initializeDSegments();
        initializeDVerts();

        // For controls, computes their corner normals. For non-controls, this method does nothing (return -1)
        int computeCornerNormals();

        // Store control points as a list
        // Store verts as first copying the curvenet. Then add middle points after that spline by spline
        // This preserves indexing for controls
        std::vector<dvert> dVerts;
        std::vector<dsegment> dSegments;
        std::vector<dspline> dSplines;
        
        
};

}   // namespace DCurvenet