// pdcurvenet.hpp
#pragma once

#include "dcurvenet.hpp"
#include "components/dvert.hpp"
#include "components/dsegment.hpp"
#include "components/dspline.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>


namespace DCurvenet {

class pdcurvenet : public dcurvenet {
    public:
        // Constructor copies in the values of the discrete curvenet, but with our modified data types
        // TODO: input mesh
        pdcurvenet(dcurvenet& parentDC);

        // TODO: Needs getters so that others can ask for internals

        // Big function for applying projection and geodesics
        // Note that eps value should be stored in mesh
        // TODO: Add mesh
        void ComputePDC();
    protected:
        // No class inheritance
    private:
        // Hard copy in values during intialization
        void copyDC(const dcurvenet& parentDC);

        // Apply projection
        void projectToSurface();

        // Snap to nearby verts/edges
        void snapToNearbyElements();

        // We copy over and then modify these lists in our PDC
        // std::vector<dvert> dVerts;
        // std::vector<dsegment> dSegments;
        // std::vector<dspline> dSplines;
        
        
};

}   // namespace DCurvenet