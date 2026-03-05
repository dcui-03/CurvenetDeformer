// cutmesh.hpp
#pragma once

#include "dcurvenet/pdcurvenet.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <map>
#include <vector>


namespace CutMesh {

class cutmesh {
    public:
        // Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
        // TODO: Add in the correct mesh data struct
        cutmesh(DCurvenet::pdcurvenet& pDC);

        // TODO: Needs getters so that the UI can ask for internals

        // Store map from cut-mesh indices to mesh indices
        std::map<int, int> cutVertToMeshVert;
    protected:
        // No class inheritance
    private:
        // "rewire" the mesh by walking across each spline around corners
        void computeCuts();

        // Mesh as a HE data structure
        // TODO: Need to pick a HE mesh class (ex. GeometryCentral or minimesh)

        // Store map from mesh vertices to their associated discrete curvenet segment(s)
        std::map<int, int> cutVertToDCSegment;
};

}   // namespace CutMesh