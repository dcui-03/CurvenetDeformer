// cutmesh.hpp
#pragma once

#include "mesh.hpp"
#include "dcurvenet/pdcurvenet.hpp"
#include "dcurvenet/components/dvert.hpp"
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

        // Recursively compute straightest geodesic between any two points
        // start is the current start point, target is the goal point,
        // Direc is the movement direction, where el_type and el_idx is what element we should walk on
        // also include a depth value which records how many recursion steps we've gone
        bool computeStraightestGeodesic(DCurvenet::dvert& start, DCurvenet::dvert& target,
                                        Eigen::Vector3d direc, int el_type, int el_idx,
                                        std::vector<DCurvenet::dvert>& dVertList, int& depth);

        // Helper function which checks start and target's adjacent elements and returns true if there is an adjacency
        bool adjacencyCheck(const DCurvenet::dvert& start, const DCurvenet::dvert& target, std::vector<int>& sharedFaces);

        // "rewire" the mesh by walking across each spline around corners
        // Don't forget to recompute element normals after! (not important, just for completeness)
        void computeCuts();

        // Inherits mesh object from parent

        // Num of each vertex that we use in the solve
        int num_V;
        int num_C;
        
};

}   // namespace Mesh