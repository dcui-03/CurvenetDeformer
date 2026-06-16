// dcurvenet.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include <Eigen/Core>

namespace DCurvenet {

// Discrete curvenet (i.e., polylines)
class dcurvenet {
    public:
        // Takes the original curvenet and discretizes it
        dcurvenet(const Curvenet::curvenet& CN);
        // Initialize with empty constructor
        dcurvenet();

    protected:
    private:
        // Reorganize during intialization
        void intializeDSplines(Curvenet::curvenet& CN);
        void initializeDSegments(Curvenet::curvenet& CN);
        void initializeDVerts(Curvenet::curvenet& CN);

        // For controls, computes their corner normals. For non-controls, this method does nothing (return -1)
        int computeCornerNormals();
        

        int sampling_param = 5;
};

}   // namespace DCurvenet