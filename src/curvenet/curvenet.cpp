#include "curvenet.hpp"

#include "components/control.hpp"
#include "components/spline.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>


namespace Curvenet {

    // Initialize from an existing list of controls and normals
    // NOTE: Initializes assuming there are already no duplicates
    curvenet::curvenet(std::vector<Eigen::Vector3d> controlP,
                       std::vector<Eigen::Vector3d> controlNormals,
                       std::vector<Eigen::Vector3d> tangentP,
                       std::vector<std::array<int, 4>> curveC) {
        // Create controls
        initializeControls(controlP, controlNormals);
        // Create tangents, then create splines
        for (int c = 0; c < curveC.size(); c++) {
            std::array<int, 2> sControl = {curveC[c][0], curveC[c][3]};
            std::array<Eigen::Vector3d, 2> sTangent = {tangentP[curveC[c][1]], tangentP[curveC[c][2]]};
            addSpline(sControl, sTangent);
        }

        // Compute the CCW spline ordering for controls
        // NOTE: no need to do this here, since we do it at spline addition
        // computeCCWOrderingAll();
    }

    // empty initializer
    curvenet::curvenet() {
        controlPoints.clear();
        tangentPoints.clear();
        splines.clear();
    }

    // Convert a curvenet object into something polyscope can read
    // Try to preserve control/tangent indexing as much as possible
    void curvenet::convertControlsToPC(std::vector<Eigen::Vector3d>& psControls) {
        psControls.clear();
        psControls.resize(controlPoints.size());
        for (int c = 0; c < controlPoints.size(); c++) {
            psControls[c] = controlPoints[c].getPos();
        }
        return;
    }
    void curvenet::convertTangentsToPC(std::vector<Eigen::Vector3d>& psTangents) {
        psTangents.clear();
        psTangents.resize(tangentPoints.size());
        for (int t = 0; t < tangentPoints.size(); t++) {
            psTangents[t] = tangentPoints[t].getPos();
        }
        return;
    }
    void curvenet::convertCurvnetToCN(std::vector<Eigen::Vector3d>& psSamples, std::vector<std::array<int, 2>>& psE) {
        psE.clear();
        // First, get all control points
        convertControlsToPC(psSamples);
        // Iterate over splines and sample them
        for (int s = 0; s < splines.size(); s++) {
            int v_idx = psSamples.size();
            std::vector<Eigen::Vector3d> UnifSamples;
            splines[s].UniformSampling(UnifSamples, false, false);
            // Start segment
            std::array<int, 2> start_segment = {splines[s].c0, v_idx};
            psE.push_back(start_segment);
            psSamples.push_back(UnifSamples[0]);
            // Interior segments
            for (int sample = 1; sample < UnifSamples.size(); sample++) {
                v_idx++;
                psSamples.push_back(UnifSamples[sample]);
                std::array<int, 2> temp_segment = {v_idx - 1, v_idx};
                psE.push_back(temp_segment);
            }
            // end segment
            std::array<int, 2> end_segment = {v_idx, splines[s].c1};
            psE.push_back(end_segment);
        }
        return;
    }
    void curvenet::convertControlsAndTangentsToCN(std::vector<Eigen::Vector3d>& psCAndT, std::vector<std::array<int, 2>>& psCAndTE) {
        psCAndTE.clear();
        // First, get all control points
        convertControlsToPC(psCAndT);
        // Iterate over splines and add their tangents
        for (int s = 0; s < splines.size(); s++) {
            int v_idx = psCAndT.size();
            // Start tangent
            psCAndT.push_back((tangentPoints[splines[s].t0]).getPos());
            std::array<int, 2> start_segment = {splines[s].c0, v_idx};
            psCAndTE.push_back(start_segment);
            // End tangent
            psCAndT.push_back((tangentPoints[splines[s].t1]).getPos());
            std::array<int, 2> end_segment = {splines[s].c1, v_idx+1};
            psCAndTE.push_back(end_segment);
        }
        return;
    }


    // Add a spline to the spline list.
    // Neither control already exists
    void curvenet::addSpline(std::array<Eigen::Vector3d, 2> controlP,
                             std::array<Eigen::Vector3d, 2> controlNormals,
                             std::array<Eigen::Vector3d, 2> tangentP) {
        int ctrl_size = controlPoints.size();
        int tan_size = tangentPoints.size();
        std::vector<Eigen::Vector3d> controlPVec(std::begin(controlP), std::end(controlP));
        std::vector<Eigen::Vector3d> normalVec(std::begin(controlNormals), std::end(controlNormals));
        initializeControls(controlPVec, normalVec);
        int spline_idx = splines.size();
        // Initialize tangents
        tangent tan0(tangentP[0], ctrl_size, spline_idx);
        tangent tan1(tangentP[1], ctrl_size + 1, spline_idx);
        tangentPoints.push_back(tan0);
        tangentPoints.push_back(tan1);
        // Initialize spline
        std::array<int, 2> ctrlIdx = {ctrl_size, ctrl_size + 1};
        std::array<int, 2> tanIdx = {tan_size, tan_size + 1};
        spline new_spline(&controlPoints, &tangentPoints, ctrlIdx, tanIdx);
        splines.push_back(new_spline);
        // Compute CCW Ordering for the 2 new controls
        controlPoints[ctrl_size].addSplineIdx(spline_idx);
        controlPoints[ctrl_size + 1].addSplineIdx(spline_idx);
        // TODO: Uncomment after adding these
        //computeCCWOrdering(ctrl_size);
        //computeCCWOrdering(ctrl_size + 1);
        return;
    }
    // Both controls already exist
    void curvenet::addSpline(std::array<int, 2> controlP, std::array<Eigen::Vector3d, 2> tangentP) {
        if (controlP[0] < 0 || controlP[0] >= static_cast<int>(controlPoints.size())) {
            throw std::out_of_range("control0 index out of range");
        }
        if (controlP[1] < 0 || controlP[1] >= static_cast<int>(controlPoints.size())) {
            throw std::out_of_range("control1 index out of range");
        }
        int tan_size = tangentPoints.size();
        int spline_idx = splines.size();
        // Initialize tangents
        tangent tan0(tangentP[0], controlP[0], spline_idx);
        tangent tan1(tangentP[1], controlP[1], spline_idx);
        tangentPoints.push_back(tan0);
        tangentPoints.push_back(tan1);
        // Initialize spline
        std::array<int, 2> tanIdx = {tan_size, tan_size + 1};
        spline new_spline(&controlPoints, &tangentPoints, controlP, tanIdx);
        splines.push_back(new_spline);
        // Compute CCW Ordering for the 2 new controls
        controlPoints[controlP[0]].addSplineIdx(spline_idx);
        controlPoints[controlP[1]].addSplineIdx(spline_idx);
        // TODO: Uncomment
        //computeCCWOrdering(controlP[0]);
        //computeCCWOrdering(controlP[1]);
        return;
    }
    // First control exists
    void curvenet::addSpline(int control0, Eigen::Vector3d control1, Eigen::Vector3d normal1, std::array<Eigen::Vector3d, 2> tangentP) {
        if (control0 < 0 || control0 >= static_cast<int>(controlPoints.size())) {
            throw std::out_of_range("control0 index out of range");
        }
        int ctrl_size = controlPoints.size();
        int tan_size = tangentPoints.size();
        addControl(control1, normal1);
        int spline_idx = splines.size();
        // Initialize tangents
        tangent tan0(tangentP[0], control0, spline_idx);
        tangent tan1(tangentP[1], ctrl_size, spline_idx);
        tangentPoints.push_back(tan0);
        tangentPoints.push_back(tan1);
        // Initialize spline
        std::array<int, 2> ctrlIdx = {control0, ctrl_size};
        std::array<int, 2> tanIdx = {tan_size, tan_size + 1};
        spline new_spline(&controlPoints, &tangentPoints, ctrlIdx, tanIdx);
        splines.push_back(new_spline);
        // Compute CCW Ordering for the 2 new controls
        controlPoints[control0].addSplineIdx(spline_idx);
        controlPoints[ctrl_size].addSplineIdx(spline_idx);
        // TODO: Uncomment
        // computeCCWOrdering(control0);
        // computeCCWOrdering(ctrl_size);
        return;
    }
    // Second control exists
    void curvenet::addSpline(Eigen::Vector3d control0, int control1, Eigen::Vector3d normal0, std::array<Eigen::Vector3d, 2> tangentP) {
        if (control1 < 0 || control1 >= static_cast<int>(controlPoints.size())) {
            throw std::out_of_range("control1 index out of range");
        }
        int ctrl_size = controlPoints.size();
        int tan_size = tangentPoints.size();
        addControl(control0, normal0);
        int spline_idx = splines.size();
        // Initialize tangents
        tangent tan0(tangentP[0], ctrl_size, spline_idx);
        tangent tan1(tangentP[1], control1, spline_idx);
        tangentPoints.push_back(tan0);
        tangentPoints.push_back(tan1);
        // Initialize spline
        std::array<int, 2> ctrlIdx = {ctrl_size, control1};
        std::array<int, 2> tanIdx = {tan_size, tan_size + 1};
        spline new_spline(&controlPoints, &tangentPoints, ctrlIdx, tanIdx);
        splines.push_back(new_spline);
        // Compute CCW Ordering for the 2 new controls
        controlPoints[ctrl_size].addSplineIdx(spline_idx);
        controlPoints[control1].addSplineIdx(spline_idx);
        // TODO: Uncomment
        //computeCCWOrdering(ctrl_size);
        //computeCCWOrdering(control1);
        return;
    }

    // Add all controls into the list
    void curvenet::initializeControls(std::vector<Eigen::Vector3d>& controlP,
                            std::vector<Eigen::Vector3d>& controlNormals) {
        for (int c = 0; c < controlP.size(); c++) {
            addControl(controlP[c], controlNormals[c]);
        }
        return;
    }

    void curvenet::addControl(Eigen::Vector3d controlP, Eigen::Vector3d controlNormal) {
        // create new control
        control new_ctrl(controlP, controlNormal);
        controlPoints.push_back(new_ctrl);
        return;
    }

    // Compute CCW Ordering of outgoing splines
    void curvenet::computeCCWOrderingAll() {
        for (int c = 0; c < controlPoints.size(); c++) {
            computeCCWOrdering(c);
        }
        return;
    }
    void curvenet::computeCCWOrdering(int ctrl_idx) {
        Eigen::Vector3d ctrl_pos = controlPoints[ctrl_idx].getPos();
        Eigen::Vector3d ctrl_Normal = controlPoints[ctrl_idx].getNormal();
        std::vector<int> adjSplines = controlPoints[ctrl_idx].getSplineIdxs();
        std::vector<Eigen::Vector2d> projectedVecs(adjSplines.size());
        std::vector<int> reordering(adjSplines.size());
        // pre-compute a basis
        Eigen::Vector3d t0, t1;
        Utils::buildPlaneBasis(ctrl_Normal, t0, t1);
        // First, get projections and write them in new basis
        for (int s = 0; s < adjSplines.size(); s++) {
            // Check which endpoint we match and grab tangent vector
            Eigen::Vector3d t;
            if (splines[adjSplines[s]].c0 == ctrl_idx) {
                t = tangentPoints[splines[adjSplines[s]].t0].getPos() - ctrl_pos;
            } else {
                t = tangentPoints[splines[adjSplines[s]].t0].getPos() - ctrl_pos;
            }
            Eigen::Vector3d proj;
            // Project onto tangent plane
            double success = Utils::projectVectorOntoTangentPlane(ctrl_Normal, t, proj);
            if (success < 0.0) {    // degenerate vector
                projectedVecs[s] = Eigen::Vector2d::Zero();
            } else {                // Nondegenerate: rewrite in new basis
                projectedVecs[s] = Utils::convertTo2D(proj, ctrl_pos, t0, t1);
            }

        }
        std::vector<int> nondegenerate;
        std::vector<int> degenerate;    // Store a list of the degenerate projections (which are arbitrary in list index)
        for (int s = 0; s < adjSplines.size(); s++) {
            // check if degenerate
            if (projectedVecs[s] == Eigen::Vector2d::Zero()) {
                degenerate.push_back(s);    // store local index
            } else {
                nondegenerate.push_back(s);
            }
        }
        // Add in the first nondegenerate vector
        if (degenerate.size() == adjSplines.size()) {
            // All projected vectors are degenerate; no need to reorder
            return;
        }
        // Iterate over the nondegenerate vertices and do signed angle
        std::vector<int> pos_increasing;    // Indices of increasing positive signed angle
        std::vector<double> pos_SA;
        std::vector<int> neg_increasing;    // Indices of increasing negative signed angle
        std::vector<double> neg_SA;
        for (int p = 1; p < nondegenerate.size(); p++) {
            // Compute signed angle compared to the first nondegenerate vector
            double dotprod = projectedVecs[nondegenerate[0]].dot(projectedVecs[nondegenerate[p]]);
            double crossprod = projectedVecs[nondegenerate[0]](0) * projectedVecs[nondegenerate[p]](1) - 
                               projectedVecs[nondegenerate[0]](1) * projectedVecs[nondegenerate[p]](0);
            double sAngle = std::atan2(crossprod, dotprod);
            // Insert into ordering based on signed angle
            if (sAngle >= 0.0) {
                std::vector<double>::iterator it = std::lower_bound(pos_SA.begin(), pos_SA.end(), sAngle);
                pos_SA.insert(it, sAngle);
                int insertion_idx = it - pos_SA.begin();
                pos_increasing.insert(pos_increasing.begin() + insertion_idx, nondegenerate[p]);
            } else {
                std::vector<double>::iterator it = std::lower_bound(neg_SA.begin(), neg_SA.end(), sAngle);
                neg_SA.insert(it, sAngle);
                int insertion_idx = it - neg_SA.begin();
                neg_increasing.insert(neg_increasing.begin() + insertion_idx, nondegenerate[p]);
            }
        }

        // Set adjacent nondegenerate splines in correct ordering
        reordering[0] = adjSplines[nondegenerate[0]];   // Fix the first one
        for (int s = 0; s < pos_SA.size(); s++) {   // Get positive signed angles (0 -> 180)
            reordering[s + 1] = adjSplines[pos_SA[s]];
        }
        for (int s = 0; s < neg_SA.size(); s++) {   // Get negative signed angles (-180 -> 0)
            reordering[pos_SA.size() + s] = adjSplines[pos_SA[s]];
        }
        for (int i = 0; i < degenerate.size(); i++) {   // Fill the rest with all the degenerate vectors
            reordering[nondegenerate.size() + i] = adjSplines[degenerate[i]];
        }

        // Set equivalence
        std::cout << ctrl_idx << " Old Spline ordering:";
        for (int s = 0; s < adjSplines.size(); s++) {
            std::cout << " " << adjSplines[s];
        }
        std::cout << std::endl;
        controlPoints[ctrl_idx].setSplineIdxs(reordering);
        std::cout << ctrl_idx << " New Spline ordering:";
        for (int s = 0; s < reordering.size(); s++) {
            std::cout << " " << reordering[s];
        }
        std::cout << std::endl;
        return;

    }

}   // namespace Curvenet