#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/point_cloud_vector_quantity.h"

#include <Eigen/Core>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "imgui.h"

#include "IO/io.hpp"
#include "curvenet/curvenet.hpp"
#include "profilemover/profilemover.hpp"
#include "utils/jsonUtils.hpp"
#include "utils/utils.hpp"

namespace {

enum class EditTarget {
    None,
    Control,
    Tangent
};

glm::vec3 toGlm(const Eigen::Vector3d& v) {
    return glm::vec3(static_cast<float>(v.x()),
                     static_cast<float>(v.y()),
                     static_cast<float>(v.z()));
}

std::vector<glm::vec3> toGlmPoints(const std::vector<Eigen::Vector3d>& points) {
    std::vector<glm::vec3> out;
    out.reserve(points.size());
    for (const auto& p : points) {
        out.push_back(toGlm(p));
    }
    return out;
}

std::vector<std::array<std::size_t, 2>> toSizeTEdges(const std::vector<std::array<int, 2>>& edges) {
    std::vector<std::array<std::size_t, 2>> out;
    out.reserve(edges.size());
    for (const auto& e : edges) {
        out.push_back({static_cast<std::size_t>(e[0]), static_cast<std::size_t>(e[1])});
    }
    return out;
}

// Mesh
std::vector<Eigen::Vector3d> meshV;
std::vector<std::vector<int>> meshF;
polyscope::SurfaceMesh* psMesh = nullptr;

// Editable curvenet data
Curvenet::curvenet neutralCurvenet;
std::unique_ptr<Curvenet::curvenet> editableCurvenet;

std::vector<glm::vec3> psEditableCurveP;
std::vector<std::array<std::size_t, 2>> psEditableCurveE;
polyscope::CurveNetwork* psEditableCurve = nullptr;

std::vector<glm::vec3> psEditableControlsP;
polyscope::PointCloud* psEditableControls = nullptr;

std::vector<glm::vec3> psEditableTangentsP;
polyscope::PointCloud* psEditableTangents = nullptr;

std::vector<glm::vec3> psEditableHandlesP;
std::vector<std::array<std::size_t, 2>> psEditableHandlesE;
polyscope::CurveNetwork* psEditableHandles = nullptr;

// Precompute cache + visuals
std::unique_ptr<ProfileMover::profilemover> PF;
bool precompValid = false;
int samplesPerMeanEdge = 5;
int uniformRefineSamples = 64;

std::vector<glm::vec3> psDCurvenetP;
std::vector<std::array<std::size_t, 2>> psDCurvenetE;
polyscope::CurveNetwork* psDCurvenet = nullptr;

std::vector<glm::vec3> psDCurvenetInteriorP;
std::vector<glm::vec3> psDCurvenetControlP;
polyscope::PointCloud* psDCurvenetInterior = nullptr;
polyscope::PointCloud* psDCurvenetControls = nullptr;

std::vector<glm::vec3> psProjectedVertexSamplesP;
std::vector<glm::vec3> psProjectedEdgeSamplesP;
std::vector<glm::vec3> psProjectedFaceSamplesP;
polyscope::PointCloud* psProjectedVertexSamples = nullptr;
polyscope::PointCloud* psProjectedEdgeSamples = nullptr;
polyscope::PointCloud* psProjectedFaceSamples = nullptr;

std::vector<glm::vec3> psControlNormalOriginsP;
std::vector<glm::vec3> psControlNormalVectors;
polyscope::PointCloud* psControlNormals = nullptr;
polyscope::PointCloudVectorQuantity* psControlNormalVectorsQ = nullptr;

std::vector<glm::vec3> psPlusSegmentFrameOriginsP;
std::vector<glm::vec3> psPlusSegmentFrameTangents;
std::vector<glm::vec3> psPlusSegmentFrameBinormals;
std::vector<glm::vec3> psPlusSegmentFrameNormals;
polyscope::PointCloud* psPlusSegmentFrames = nullptr;
polyscope::PointCloudVectorQuantity* psPlusSegmentFrameTangentsQ = nullptr;
polyscope::PointCloudVectorQuantity* psPlusSegmentFrameBinormalsQ = nullptr;
polyscope::PointCloudVectorQuantity* psPlusSegmentFrameNormalsQ = nullptr;

std::vector<glm::vec3> psMinusSegmentFrameOriginsP;
std::vector<glm::vec3> psMinusSegmentFrameTangents;
std::vector<glm::vec3> psMinusSegmentFrameBinormals;
std::vector<glm::vec3> psMinusSegmentFrameNormals;
polyscope::PointCloud* psMinusSegmentFrames = nullptr;
polyscope::PointCloudVectorQuantity* psMinusSegmentFrameTangentsQ = nullptr;
polyscope::PointCloudVectorQuantity* psMinusSegmentFrameBinormalsQ = nullptr;
polyscope::PointCloudVectorQuantity* psMinusSegmentFrameNormalsQ = nullptr;

// UI state
bool editControlsMode = false;
bool editTangentsMode = false;
bool constrainTangentsToPlane = true;
bool showCornerNormals = false;
bool showSegmentFrames = false;
EditTarget activeTarget = EditTarget::None;
int activeIndex = -1;
bool gizmoActive = false;
polyscope::TransformationGizmo* vertexGizmo = nullptr;

std::string inputPath;
std::string curvesPath;

void updatePrecomputeVisibility() {
    const bool enabled = precompValid;

    if (psDCurvenet != nullptr) psDCurvenet->setEnabled(enabled);
    if (psDCurvenetInterior != nullptr) psDCurvenetInterior->setEnabled(enabled);
    if (psDCurvenetControls != nullptr) psDCurvenetControls->setEnabled(enabled);
    if (psProjectedVertexSamples != nullptr) psProjectedVertexSamples->setEnabled(enabled);
    if (psProjectedEdgeSamples != nullptr) psProjectedEdgeSamples->setEnabled(enabled);
    if (psProjectedFaceSamples != nullptr) psProjectedFaceSamples->setEnabled(enabled);

    const bool normalsEnabled = enabled && showCornerNormals;
    if (psControlNormals != nullptr) psControlNormals->setEnabled(normalsEnabled);
    if (psControlNormalVectorsQ != nullptr) psControlNormalVectorsQ->setEnabled(normalsEnabled);

    const bool framesEnabled = enabled && showSegmentFrames;
    if (psPlusSegmentFrames != nullptr) psPlusSegmentFrames->setEnabled(framesEnabled);
    if (psPlusSegmentFrameTangentsQ != nullptr) psPlusSegmentFrameTangentsQ->setEnabled(framesEnabled);
    if (psPlusSegmentFrameBinormalsQ != nullptr) psPlusSegmentFrameBinormalsQ->setEnabled(framesEnabled);
    if (psPlusSegmentFrameNormalsQ != nullptr) psPlusSegmentFrameNormalsQ->setEnabled(framesEnabled);
    if (psMinusSegmentFrames != nullptr) psMinusSegmentFrames->setEnabled(framesEnabled);
    if (psMinusSegmentFrameTangentsQ != nullptr) psMinusSegmentFrameTangentsQ->setEnabled(framesEnabled);
    if (psMinusSegmentFrameBinormalsQ != nullptr) psMinusSegmentFrameBinormalsQ->setEnabled(framesEnabled);
    if (psMinusSegmentFrameNormalsQ != nullptr) psMinusSegmentFrameNormalsQ->setEnabled(framesEnabled);
}

void removeGizmo() {
    if (vertexGizmo != nullptr) {
        vertexGizmo->remove();
        vertexGizmo = nullptr;
    }
    gizmoActive = false;
}

void clearActiveEdit() {
    activeTarget = EditTarget::None;
    activeIndex = -1;
    removeGizmo();
}

void addGizmoAt(const Eigen::Vector3d& position) {
    gizmoActive = true;
    if (vertexGizmo == nullptr) {
        vertexGizmo = polyscope::addTransformationGizmo("curve_edit_gizmo");
        vertexGizmo->setAllowTranslation(true);
        vertexGizmo->setAllowRotation(false);
        vertexGizmo->setAllowScaling(false);
        vertexGizmo->setInteractInLocalSpace(false);
    }
    vertexGizmo->setPosition(toGlm(position));
}

void hidePrecomputeVisuals() {
    updatePrecomputeVisibility();
}

void invalidatePrecompute() {
    precompValid = false;
    PF.reset();
    hidePrecomputeVisuals();
}

void refreshEditableCurvenetVisuals() {
    if (!editableCurvenet) {
        return;
    }

    std::vector<Eigen::Vector3d> curveSamples;
    std::vector<std::array<int, 2>> curveEdges;
    editableCurvenet->convertCurvnetToCN(curveSamples, curveEdges, 96);
    psEditableCurveP = toGlmPoints(curveSamples);
    psEditableCurveE = toSizeTEdges(curveEdges);
    psEditableCurve = polyscope::registerCurveNetwork("Editable Curvenet", psEditableCurveP, psEditableCurveE);
    psEditableCurve->setColor({0.92f, 0.55f, 0.14f});
    psEditableCurve->setRadius(0.0025f, true);
    psEditableCurve->setEnabled(true);

    std::vector<Eigen::Vector3d> controls;
    editableCurvenet->convertControlsToPC(controls);
    psEditableControlsP = toGlmPoints(controls);
    psEditableControls = polyscope::registerPointCloud("Editable Controls", psEditableControlsP);
    psEditableControls->setPointColor({0.95f, 0.20f, 0.12f});
    psEditableControls->setPointRadius(0.0080, true);
    psEditableControls->setEnabled(true);

    std::vector<Eigen::Vector3d> tangents;
    editableCurvenet->convertTangentsToPC(tangents);
    psEditableTangentsP = toGlmPoints(tangents);
    psEditableTangents = polyscope::registerPointCloud("Editable Tangents", psEditableTangentsP);
    psEditableTangents->setPointColor({0.12f, 0.82f, 0.30f});
    psEditableTangents->setPointRadius(0.0060, true);
    psEditableTangents->setEnabled(true);

    std::vector<Eigen::Vector3d> handles;
    std::vector<std::array<int, 2>> handleEdges;
    editableCurvenet->convertControlsAndTangentsToCN(handles, handleEdges);
    psEditableHandlesP = toGlmPoints(handles);
    psEditableHandlesE = toSizeTEdges(handleEdges);
    psEditableHandles = polyscope::registerCurveNetwork("Editable Handles", psEditableHandlesP, psEditableHandlesE);
    psEditableHandles->setColor({0.55f, 0.70f, 0.12f});
    psEditableHandles->setRadius(0.0012f, true);
    psEditableHandles->setEnabled(true);
}

void refreshPrecomputeVisuals() {
    if (!PF) {
        hidePrecomputeVisuals();
        return;
    }

    std::vector<std::array<double, 3>> dcurvePtsD;
    std::vector<std::array<std::size_t, 2>> dcurveEdgesD;
    Utils::buildPolyscopeDiscreteCurveNetwork(PF->getNeutralDCurvenet(), dcurvePtsD, dcurveEdgesD);
    psDCurvenetP.clear();
    psDCurvenetP.reserve(dcurvePtsD.size());
    for (const auto& p : dcurvePtsD) {
        psDCurvenetP.emplace_back(static_cast<float>(p[0]),
                                  static_cast<float>(p[1]),
                                  static_cast<float>(p[2]));
    }
    psDCurvenetE = dcurveEdgesD;
    psDCurvenet = polyscope::registerCurveNetwork("Precompute: DCurveNet", psDCurvenetP, psDCurvenetE);
    psDCurvenet->setColor({0.15f, 0.80f, 0.95f});
    psDCurvenet->setRadius(0.0040f, true);
    psDCurvenet->setEnabled(true);

    psDCurvenetInteriorP.clear();
    psDCurvenetControlP.clear();
    const auto& dVerts = PF->getNeutralDCurvenet().verts();
    for (const auto& v : dVerts) {
        const glm::vec3 gp = toGlm(v.position());
        if (v.isControl()) {
            psDCurvenetControlP.push_back(gp);
        } else {
            psDCurvenetInteriorP.push_back(gp);
        }
    }
    psDCurvenetInterior = polyscope::registerPointCloud("Precompute: DCurveNet Interior", psDCurvenetInteriorP);
    psDCurvenetInterior->setPointColor({0.60f, 0.92f, 1.00f});
    psDCurvenetInterior->setPointRadius(0.0035, true);
    psDCurvenetInterior->setEnabled(true);

    psDCurvenetControls = polyscope::registerPointCloud("Precompute: DCurveNet Controls", psDCurvenetControlP);
    psDCurvenetControls->setPointColor({1.00f, 0.42f, 0.10f});
    psDCurvenetControls->setPointRadius(0.0065, true);
    psDCurvenetControls->setEnabled(true);

    psProjectedVertexSamplesP.clear();
    psProjectedEdgeSamplesP.clear();
    psProjectedFaceSamplesP.clear();
    for (const auto& sample : PF->getNeutralPDCurvenet().projectedSamples()) {
        const glm::vec3 gp = toGlm(sample.projected_position);
        switch (sample.attachment) {
            case DCurvenet::MeshAttachmentType::Vertex:
                psProjectedVertexSamplesP.push_back(gp);
                break;
            case DCurvenet::MeshAttachmentType::Edge:
                psProjectedEdgeSamplesP.push_back(gp);
                break;
            case DCurvenet::MeshAttachmentType::Face:
                psProjectedFaceSamplesP.push_back(gp);
                break;
        }
    }

    psProjectedVertexSamples = polyscope::registerPointCloud("Precompute: Projected Vertex Samples", psProjectedVertexSamplesP);
    psProjectedVertexSamples->setPointColor({1.00f, 0.25f, 0.20f});
    psProjectedVertexSamples->setPointRadius(0.0090, true);
    psProjectedVertexSamples->setEnabled(true);

    psProjectedEdgeSamples = polyscope::registerPointCloud("Precompute: Projected Edge Samples", psProjectedEdgeSamplesP);
    psProjectedEdgeSamples->setPointColor({1.00f, 0.85f, 0.20f});
    psProjectedEdgeSamples->setPointRadius(0.0080, true);
    psProjectedEdgeSamples->setEnabled(true);

    psProjectedFaceSamples = polyscope::registerPointCloud("Precompute: Projected Face Samples", psProjectedFaceSamplesP);
    psProjectedFaceSamples->setPointColor({0.15f, 0.95f, 0.45f});
    psProjectedFaceSamples->setPointRadius(0.0070, true);
    psProjectedFaceSamples->setEnabled(true);

    std::vector<std::array<double, 3>> controlNormalOriginsD;
    std::vector<std::array<double, 3>> controlNormalVectorsD;
    Utils::buildPolyscopeControlCornerNormals(PF->getNeutralDCurvenet(),
                                              controlNormalOriginsD,
                                              controlNormalVectorsD);
    psControlNormalOriginsP.clear();
    psControlNormalVectors.clear();
    for (const auto& p : controlNormalOriginsD) {
        psControlNormalOriginsP.emplace_back(static_cast<float>(p[0]),
                                             static_cast<float>(p[1]),
                                             static_cast<float>(p[2]));
    }
    for (const auto& v : controlNormalVectorsD) {
        psControlNormalVectors.emplace_back(static_cast<float>(v[0]),
                                            static_cast<float>(v[1]),
                                            static_cast<float>(v[2]));
    }
    psControlNormals = polyscope::registerPointCloud("Precompute: Corner Normal Origins", psControlNormalOriginsP);
    psControlNormals->setPointColor({0.10f, 0.95f, 0.25f});
    psControlNormals->setPointRadius(0.0050, true);
    psControlNormals->setEnabled(false);
    psControlNormalVectorsQ = nullptr;
    if (!psControlNormalVectors.empty()) {
        psControlNormalVectorsQ = psControlNormals->addVectorQuantity("Corner Normals",
                                                                     psControlNormalVectors,
                                                                     polyscope::VectorType::STANDARD);
        psControlNormalVectorsQ->setVectorColor({0.10f, 0.95f, 0.25f});
        psControlNormalVectorsQ->setVectorLengthRange(1.0);
        psControlNormalVectorsQ->setVectorLengthScale(0.08, false);
        psControlNormalVectorsQ->setEnabled(false);
    }

    psPlusSegmentFrameOriginsP.clear();
    psPlusSegmentFrameTangents.clear();
    psPlusSegmentFrameBinormals.clear();
    psPlusSegmentFrameNormals.clear();
    psMinusSegmentFrameOriginsP.clear();
    psMinusSegmentFrameTangents.clear();
    psMinusSegmentFrameBinormals.clear();
    psMinusSegmentFrameNormals.clear();
    const auto& segs = PF->getNeutralDCurvenet().segments();
    for (const auto& seg : segs) {
        const Eigen::Vector3d midpoint =
            0.5 * (dVerts[static_cast<std::size_t>(seg.startDvert())].position() +
                   dVerts[static_cast<std::size_t>(seg.endDvert())].position());
        if (seg.plusSide().valid) {
            psPlusSegmentFrameOriginsP.push_back(toGlm(midpoint));
            psPlusSegmentFrameTangents.push_back(toGlm(seg.plusSide().BS.col(0)));
            psPlusSegmentFrameBinormals.push_back(toGlm(seg.plusSide().BS.col(1)));
            psPlusSegmentFrameNormals.push_back(toGlm(seg.plusSide().BS.col(2)));
        }
        if (seg.minusSide().valid) {
            psMinusSegmentFrameOriginsP.push_back(toGlm(midpoint));
            psMinusSegmentFrameTangents.push_back(toGlm(seg.minusSide().BS.col(0)));
            psMinusSegmentFrameBinormals.push_back(toGlm(seg.minusSide().BS.col(1)));
            psMinusSegmentFrameNormals.push_back(toGlm(seg.minusSide().BS.col(2)));
        }
    }

    psPlusSegmentFrames = polyscope::registerPointCloud("Precompute: + Frame Origins", psPlusSegmentFrameOriginsP);
    psPlusSegmentFrames->setPointColor({0.92f, 0.92f, 0.92f});
    psPlusSegmentFrames->setPointRadius(0.0035, true);
    psPlusSegmentFrames->setEnabled(false);
    psPlusSegmentFrameTangentsQ = nullptr;
    psPlusSegmentFrameBinormalsQ = nullptr;
    psPlusSegmentFrameNormalsQ = nullptr;
    if (!psPlusSegmentFrameTangents.empty()) {
        psPlusSegmentFrameTangentsQ = psPlusSegmentFrames->addVectorQuantity("+ Frame Tangent",
                                                                            psPlusSegmentFrameTangents,
                                                                            polyscope::VectorType::STANDARD);
        psPlusSegmentFrameTangentsQ->setVectorColor({0.95f, 0.25f, 0.20f});
        psPlusSegmentFrameTangentsQ->setVectorLengthRange(1.0);
        psPlusSegmentFrameTangentsQ->setVectorLengthScale(1.0, false);
        psPlusSegmentFrameTangentsQ->setEnabled(false);

        psPlusSegmentFrameBinormalsQ = psPlusSegmentFrames->addVectorQuantity("+ Frame Binormal",
                                                                             psPlusSegmentFrameBinormals,
                                                                             polyscope::VectorType::STANDARD);
        psPlusSegmentFrameBinormalsQ->setVectorColor({0.98f, 0.75f, 0.20f});
        psPlusSegmentFrameBinormalsQ->setVectorLengthRange(1.0);
        psPlusSegmentFrameBinormalsQ->setVectorLengthScale(1.0, false);
        psPlusSegmentFrameBinormalsQ->setEnabled(false);

        psPlusSegmentFrameNormalsQ = psPlusSegmentFrames->addVectorQuantity("+ Frame Normal",
                                                                           psPlusSegmentFrameNormals,
                                                                           polyscope::VectorType::STANDARD);
        psPlusSegmentFrameNormalsQ->setVectorColor({0.12f, 0.95f, 0.35f});
        psPlusSegmentFrameNormalsQ->setVectorLengthRange(1.0);
        psPlusSegmentFrameNormalsQ->setVectorLengthScale(1.0, false);
        psPlusSegmentFrameNormalsQ->setEnabled(false);
    }

    psMinusSegmentFrames = polyscope::registerPointCloud("Precompute: - Frame Origins", psMinusSegmentFrameOriginsP);
    psMinusSegmentFrames->setPointColor({0.75f, 0.75f, 0.75f});
    psMinusSegmentFrames->setPointRadius(0.0035, true);
    psMinusSegmentFrames->setEnabled(false);
    psMinusSegmentFrameTangentsQ = nullptr;
    psMinusSegmentFrameBinormalsQ = nullptr;
    psMinusSegmentFrameNormalsQ = nullptr;
    if (!psMinusSegmentFrameTangents.empty()) {
        psMinusSegmentFrameTangentsQ = psMinusSegmentFrames->addVectorQuantity("- Frame Tangent",
                                                                              psMinusSegmentFrameTangents,
                                                                              polyscope::VectorType::STANDARD);
        psMinusSegmentFrameTangentsQ->setVectorColor({0.85f, 0.45f, 0.65f});
        psMinusSegmentFrameTangentsQ->setVectorLengthRange(1.0);
        psMinusSegmentFrameTangentsQ->setVectorLengthScale(1.0, false);
        psMinusSegmentFrameTangentsQ->setEnabled(false);

        psMinusSegmentFrameBinormalsQ = psMinusSegmentFrames->addVectorQuantity("- Frame Binormal",
                                                                               psMinusSegmentFrameBinormals,
                                                                               polyscope::VectorType::STANDARD);
        psMinusSegmentFrameBinormalsQ->setVectorColor({0.45f, 0.80f, 0.95f});
        psMinusSegmentFrameBinormalsQ->setVectorLengthRange(1.0);
        psMinusSegmentFrameBinormalsQ->setVectorLengthScale(1.0, false);
        psMinusSegmentFrameBinormalsQ->setEnabled(false);

        psMinusSegmentFrameNormalsQ = psMinusSegmentFrames->addVectorQuantity("- Frame Normal",
                                                                             psMinusSegmentFrameNormals,
                                                                             polyscope::VectorType::STANDARD);
        psMinusSegmentFrameNormalsQ->setVectorColor({0.35f, 0.65f, 1.00f});
        psMinusSegmentFrameNormalsQ->setVectorLengthRange(1.0);
        psMinusSegmentFrameNormalsQ->setVectorLengthScale(1.0, false);
        psMinusSegmentFrameNormalsQ->setEnabled(false);
    }

    updatePrecomputeVisibility();
}

void resetEditableCurvenet() {
    if (!editableCurvenet) {
        return;
    }
    *editableCurvenet = neutralCurvenet;
    clearActiveEdit();
    refreshEditableCurvenetVisuals();
    invalidatePrecompute();
}

void runPrecompute() {
    if (!editableCurvenet) {
        return;
    }
    PF = std::make_unique<ProfileMover::profilemover>(meshV, meshF, *editableCurvenet);
    PF->setDiscretizationParameters(samplesPerMeanEdge, uniformRefineSamples);
    PF->precomputation();
    precompValid = true;
    refreshPrecomputeVisuals();
}

void markCurveEdited() {
    refreshEditableCurvenetVisuals();
    invalidatePrecompute();
}

void applyGizmoEdit() {
    if (!gizmoActive || vertexGizmo == nullptr || !editableCurvenet) {
        return;
    }

    Eigen::Vector3d newPos(vertexGizmo->getPosition().x,
                           vertexGizmo->getPosition().y,
                           vertexGizmo->getPosition().z);

    constexpr double kMoveEps = 1e-8;

    if (activeTarget == EditTarget::Tangent && activeIndex >= 0) {
        const Curvenet::tangent& t = editableCurvenet->tangentAt(activeIndex);
        const int parentControl = t.getParentControl();
        if (constrainTangentsToPlane && parentControl >= 0) {
            const Curvenet::control& c = editableCurvenet->controlAt(parentControl);
            newPos = Utils::projectPointOntoPlane(c.getNormal(), c.getPosition(), newPos);
            vertexGizmo->setPosition(toGlm(newPos));
        }

        if ((newPos - t.getPosition()).norm() > kMoveEps) {
            editableCurvenet->setTangentPosition(activeIndex, newPos);
            markCurveEdited();
        }
        return;
    }

    if (activeTarget == EditTarget::Control && activeIndex >= 0) {
        auto& controls = editableCurvenet->controlsMutable();
        auto& tangents = editableCurvenet->tangentsMutable();
        Eigen::Vector3d oldPos = controls[static_cast<std::size_t>(activeIndex)].getPosition();
        Eigen::Vector3d delta = newPos - oldPos;
        if (delta.norm() > kMoveEps) {
            controls[static_cast<std::size_t>(activeIndex)].setPosition(newPos);
            for (auto& t : tangents) {
                if (t.getParentControl() == activeIndex) {
                    t.setPosition(t.getPosition() + delta);
                }
            }
            editableCurvenet->recomputeAllSplines();
            markCurveEdited();
        }
    }
}

void myCallback() {
    const bool mouseClicked = ImGui::IsMouseClicked(0);
    const glm::vec2 screen{ImGui::GetIO().MousePos.x, ImGui::GetIO().MousePos.y};
    const polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen);

    ImGui::SeparatorText("Curve Edit");
    if (ImGui::Button(editControlsMode ? "Stop Editing Controls" : "Edit Controls")) {
        editControlsMode = !editControlsMode;
        editTangentsMode = false;
        clearActiveEdit();
    }
    if (ImGui::Button(editTangentsMode ? "Stop Editing Tangents" : "Edit Tangents")) {
        editTangentsMode = !editTangentsMode;
        editControlsMode = false;
        clearActiveEdit();
    }
    ImGui::Checkbox("Constrain Tangents To Plane", &constrainTangentsToPlane);

    if (ImGui::Button("Reset Curve")) {
        resetEditableCurvenet();
    }

    ImGui::SeparatorText("Precompute");
    if (ImGui::SliderInt("Samples / Mean Edge", &samplesPerMeanEdge, 1, 32) && precompValid) {
        runPrecompute();
    }
    if (ImGui::SliderInt("Uniform Refine Samples", &uniformRefineSamples, 8, 512) && precompValid) {
        runPrecompute();
    }

    if (ImGui::Button("Run Precompute")) {
        runPrecompute();
    }
    if (ImGui::Button(showCornerNormals ? "Hide Corner Normals" : "Show Corner Normals")) {
        showCornerNormals = !showCornerNormals;
        updatePrecomputeVisibility();
    }
    if (ImGui::Button(showSegmentFrames ? "Hide Segment Frames" : "Show Segment Frames")) {
        showSegmentFrames = !showSegmentFrames;
        updatePrecomputeVisibility();
    }

    if (ImGui::Button("Compute Deformation")) {
        std::cout << "Deformation mode is not implemented yet." << std::endl;
    }

    if (editableCurvenet) {
        ImGui::Text("Controls: %zu", editableCurvenet->controls().size());
        ImGui::Text("Tangents: %zu", editableCurvenet->tangents().size());
        ImGui::Text("Splines: %zu", editableCurvenet->getSplines().size());
    }
    if (PF && precompValid) {
        ImGui::Text("dVerts: %zu", PF->getNeutralDCurvenet().verts().size());
        ImGui::Text("dSegments: %zu", PF->getNeutralDCurvenet().segments().size());
    }

    if (mouseClicked && editControlsMode && pick.isHit && pick.structure == psEditableControls) {
        const polyscope::PointCloudPickResult pcPick = psEditableControls->interpretPickResult(pick);
        activeTarget = EditTarget::Control;
        activeIndex = static_cast<int>(pcPick.index);
        addGizmoAt(editableCurvenet->controlAt(activeIndex).getPosition());
    } else if (mouseClicked && editTangentsMode && pick.isHit && pick.structure == psEditableTangents) {
        const polyscope::PointCloudPickResult pcPick = psEditableTangents->interpretPickResult(pick);
        activeTarget = EditTarget::Tangent;
        activeIndex = static_cast<int>(pcPick.index);
        addGizmoAt(editableCurvenet->tangentAt(activeIndex).getPosition());
    } else if (mouseClicked && (editControlsMode || editTangentsMode) &&
               (!pick.isHit || (pick.structure != psEditableControls && pick.structure != psEditableTangents))) {
        clearActiveEdit();
    }

    applyGizmoEdit();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: ./profile_mover <mesh.obj> [curves.json]" << std::endl;
        return 1;
    }

    inputPath = argv[1];
    curvesPath = (argc > 2) ? argv[2] : "data/sphere-curves.json";

    if (!IO::readOBJ(inputPath, meshV, meshF)) {
        std::cerr << "Failed to open OBJ: " << inputPath << std::endl;
        return 1;
    }
    const JSONUtils::CurvenetInput cnInput = JSONUtils::loadBezierCurvenetInput(curvesPath);

    neutralCurvenet = Curvenet::curvenet(cnInput.controlP, cnInput.surfaceN, cnInput.curveC);
    editableCurvenet = std::make_unique<Curvenet::curvenet>(neutralCurvenet);

    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::state::userCallback = myCallback;
    polyscope::init();

    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", meshV, meshF);
    psMesh->setSurfaceColor({0.80f, 0.84f, 0.88f});

    refreshEditableCurvenetVisuals();
    hidePrecomputeVisuals();

    polyscope::show();
    return 0;
}
