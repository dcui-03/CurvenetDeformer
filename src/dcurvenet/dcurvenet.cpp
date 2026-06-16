#include "dcurvenet.hpp"

#include <Eigen/Core>
#include <vector>

namespace DCurvenet {

    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(const Curvenet::curvenet& CN, int alpha) {
        // 1. First copy in the control points
        //    a. Also copy over the size of the adjHE such that we can get a 1:1 match on the outgoing spline indices

        // 2. Initialize new dCN copies of the CN splines

        // 3. Iterate over each CN spline.
        //    a. Compute the discrete sampling
        //    b. Initialize new Vert, HalfEdges, Edges
        //    c. Re-wire as necessary

        // 4. Compute all corner normals
        // 5. Transport normals along all splines
        // 6. Compute scaled frames on all splines
    }

}   // namespace DCurvenet