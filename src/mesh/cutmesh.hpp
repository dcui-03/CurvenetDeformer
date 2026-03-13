// cutmesh.hpp
#pragma once

#include "mesh.hpp"
#include "dcurvenet/pdcurvenet.hpp"
#include <Eigen/Core>
#include <map>
#include <vector>


namespace Mesh {

class cutmesh : public mesh {
    public:
        // Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
        // TODO: Add in the correct mesh data struct
        cutmesh(DCurvenet::pdcurvenet& pDC);

        // TODO: Needs getters so that the UI can ask for internals

        // Store map from cut-mesh indices to mesh indices
        std::map<int, int> cutVertToMeshVert;
        // Store map from mesh vertices to their associated discrete curvenet segment(s)
        std::map<int, int> cutVertToDCSegment;
    protected:
        // No class inheritance
    private:
        // Copy in elements from input mesh
        void copyFromMesh();

        // "rewire" the mesh by walking across each spline around corners
        // Don't forget to recompute element normals after! (not important, just for completeness)
        void computeCuts();

        // Inherits mesh object from parent

        // Num of each vertex that we use in the solve
        int num_V;
        int num_C;
        
};

}   // namespace Mesh