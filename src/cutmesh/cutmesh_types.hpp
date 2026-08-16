// cutmesh_types.hpp
#pragma once

#include "mesh/mesh_types.hpp"
#include <Eigen/Core>

// File with basic structs used by cutmesh class

namespace Mesh {

    // Data needed if we want to interpolate deformations
    struct vertDeformData {
        // Vector from the projected point on the rest mesh to the rest curvenet vert (if label == 1)
        Eigen::Vector3d projVector = Eigen::Vector3d::Zero();
        // Deformation gradient eventually computed using Laplacian
        Eigen::Matrix3d defGrad = Eigen::Matrix3d::Identity();
    };

    // Extra per-vertex data only cut-vertices need, parallel to cutmesh's V (same size, same indexing)
    struct CutData {
        int label = 0;  // {0 if original mesh vertex, 1 if projected CN vertex, 2 otherwise}
        int corner_idx = -1;   // Corresponding dCN HALFEDGE index for cut-mesh (if cut-vertex is associated with a dCN vert or HE)

        vertProjData projData;  // Where this cut-vertex maps to on the reference mesh
        vertDeformData defData;
    };
}   // namespace Mesh
