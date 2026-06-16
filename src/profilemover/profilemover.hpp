// profilemover.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include "mesh/cutmesh.hpp"
#include "utils/decUtils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace ProfileMover {

class profilemover {
    public:
        // Constructor, which first builds the mesh
        profilemover(std::vector<Eigen::Vector3d>& meshV, std::vector<std::vector<int>>& meshT);

        // Precompute cut-mesh and operators
        void precomputation();

        // Apply deformation given the new control point locations (connectivity should be same)
        void deformation(std::vector<Eigen::Vector3d>& controlP);

        // TODO: Initialize directly from existing mesh data struct (ex. GeometryCentral or minimesh)

        // TODO: Needs getters so that the UI can ask for internals
    protected:
        // No class inheritance
    private:
        // initialization functions called by constructor
        // TODO: We should probably compute vertex, edge, and face normals here and store them somewhere...
        void initializeMesh(std::vector<Eigen::Vector3d>& meshV, std::vector<std::vector<int>>& meshT);
        void initializeCurvenet(std::vector<Eigen::Vector3d> controlP, std::vector<std::vector<int>> curveC);

        // TODO: Figure out what to do with the num samples func
        int spline::computeNumSamples(int alpha, double meanE, double arclength) {
            return std::max(2, static_cast<int>(alpha * (arclength)/meanE));
        }


        // Computes deformations on the temporary discrete curvenet, then fills in a list of flattened deformation
        // gradients per vertex
        void computeDeformationGradients(DCurvenet::dcurvenet& tempDC, std::vector<Eigen::VectorXd>& defGrads);

        // Store copy of initial mesh
        Mesh::mesh cutMesh;
        
        

        // NOTE: no need to store updated states curvenet and dcurvenet, we have to create new copies at execution time
        // Store the neutral curvenet (spline)
        // TODO: Do we need the original curvenet? --> Only for resets
        // We do need the original nDCurvenet so that we can compute deformation gradients
        Curvenet::curvenet nCurvenet;
        // Store the neutral discrete curvenet
        DCurvenet::dcurvenet nDCurvenet;
        // Store neutral cut-mesh
        //CutMesh::cutmesh nCutmesh;

        
        // Store operators
        // TODO: Need functions to compute V and C
        // TODO: Instead of storing VtLV, store its factorization.
        Eigen::SparseMatrix<double> VtLV;
        Eigen::SparseMatrix<double> mVtL;
        Eigen::SparseMatrix<double> C;
        
};

}   // namespace ProfileMover