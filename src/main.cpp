#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <tuple>
#include <array>
#include <vector>
#include <stdexcept>
#include <memory>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

//#include "args/args.hxx"
#include "imgui.h"

// My files
#include "profileformer/profileformer.hpp"
#include "utils/jsonUtils.hpp"
#include "utils/utils.hpp"

// Main file for visualization with Polyscope

/*
NOTES:
- We store two copies of each object. One is for Polyscope to mess with, the other is the "neutral"/original copy
*/

// VARIABLES FOR POLYSCOPE OBJECTS
// Handles for the surface mesh
std::vector<Eigen::Vector3d> psV; // Vertex list
std::vector<std::vector<int>> psT; // Face list: Note the inner list has arbitrary size for non-triangle faces
polyscope::SurfaceMesh* psMesh;

// Handles for curvenet and discrete curvenet
std::vector<glm::vec3> psCurvenetP; // Point list (smooth curvenet)
std::vector<std::array<size_t, 2>> psCurvenetE; // Edge list (smooth curvenet)
polyscope::CurveNetwork* psCurvenet = nullptr;

std::vector<glm::vec3> psDCurvenetP; // Point list (discrete curvenet)
std::vector<std::array<size_t, 2>> psDCurvenetE; // Edge list (discrete curvenet)
polyscope::CurveNetwork* psDCurvenet = nullptr;

std::vector<glm::vec3> psDCurvenetInteriorP; // Interior dcurvenet vertices
std::vector<glm::vec3> psDCurvenetControlP; // Control dcurvenet vertices
polyscope::PointCloud* psDCurvenetInterior = nullptr;
polyscope::PointCloud* psDCurvenetControls = nullptr;

std::vector<glm::vec3> psProjectedVertexSamplesP;
std::vector<glm::vec3> psProjectedEdgeSamplesP;
std::vector<glm::vec3> psProjectedFaceSamplesP;
polyscope::PointCloud* psProjectedVertexSamples = nullptr;
polyscope::PointCloud* psProjectedEdgeSamples = nullptr;
polyscope::PointCloud* psProjectedFaceSamples = nullptr;

std::vector<glm::vec3> psControlNormalOriginsP; // Point list (origins for corner normals)
std::vector<glm::vec3> psControlNormalVectors; // Vector list (corner normals)
polyscope::PointCloud* psControlNormals = nullptr;

// Handles for controls (?)
// Not sure what we need visualization-wise
std::vector<glm::vec3> psControlsP; // Point list
polyscope::PointCloud* psControls = nullptr;


// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;
std::string CurvesPath;

// Core pipeline object
std::unique_ptr<ProfileFormer::profileformer> PF;

// UI HELPERS
int selectedVertex = -1;    // Index of selected vertex on mesh
bool controlMode = false;  // Allows users to place control points
bool gizmoMode = false; // This allows the user to create a gizmo

// Gizmo helpers
bool activeGizmo = false; // This tells us if there is an active gizmo
static polyscope::TransformationGizmo* vertexGizmo = nullptr;

// INTERNAL OBJECT COPIES
// Internal mesh copy (Eigen) --> May want to use GeometryCentral or Minimesh
// Neutral mesh state
std::vector<Eigen::Vector3d> V;
std::vector<std::vector<int>> T;

// Neutral discrete curvenet
std::vector<Eigen::Vector3d> curvenetP;
std::vector<std::vector<int>> curvenetE;

// Neutral pointcloud
std::vector<Eigen::Vector3d> controlsP;
std::vector<std::vector<int>> controlsE;

// Netural Cut-mesh and Projected DC will also go here once we have them

// Forward declarations for legacy callback aliases
void removeGizmo();
void resetMeshVertexPositions();


// ----------------- FUNCTIONS BEGIN HERE -------------------------

void refreshCurvenetVisuals() {
    if (!PF) return;
    auto toGlm = [](const Eigen::Vector3d& v) {
        return glm::vec3(static_cast<float>(v.x()),
                         static_cast<float>(v.y()),
                         static_cast<float>(v.z()));
    };

    std::vector<std::array<double, 3>> curvePtsD;
    std::vector<std::array<std::size_t, 2>> curveEdgesD;
    Utils::buildPolyscopeCurveNetwork(PF->getNeutralCurvenet(), curvePtsD, curveEdgesD, 180);
    psCurvenetP.clear();
    psCurvenetP.reserve(curvePtsD.size());
    for (const auto& p : curvePtsD) {
        psCurvenetP.emplace_back(static_cast<float>(p[0]),
                                 static_cast<float>(p[1]),
                                 static_cast<float>(p[2]));
    }
    psCurvenetE = curveEdgesD;

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

    psDCurvenetInteriorP.clear();
    psDCurvenetControlP.clear();
    const auto& dVerts = PF->getNeutralDCurvenet().verts();
    psDCurvenetInteriorP.reserve(dVerts.size());
    psDCurvenetControlP.reserve(dVerts.size());
    for (const auto& v : dVerts) {
        const auto& p = v.position();
        const glm::vec3 gp(static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()));
        if (v.isControl()) {
            psDCurvenetControlP.push_back(gp);
        } else {
            psDCurvenetInteriorP.push_back(gp);
        }
    }

    psProjectedVertexSamplesP.clear();
    psProjectedEdgeSamplesP.clear();
    psProjectedFaceSamplesP.clear();
    const auto& projectedSamples = PF->getNeutralPDCurvenet().projectedSamples();
    psProjectedVertexSamplesP.reserve(projectedSamples.size());
    psProjectedEdgeSamplesP.reserve(projectedSamples.size());
    psProjectedFaceSamplesP.reserve(projectedSamples.size());
    for (const auto& sample : projectedSamples) {
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

    std::vector<std::array<double, 3>> controlNormalOriginsD;
    std::vector<std::array<double, 3>> controlNormalVectorsD;
    Utils::buildPolyscopeControlCornerNormals(PF->getNeutralDCurvenet(),
                                              controlNormalOriginsD,
                                              controlNormalVectorsD);

    psControlNormalOriginsP.clear();
    psControlNormalVectors.clear();
    psControlNormalOriginsP.reserve(controlNormalOriginsD.size());
    psControlNormalVectors.reserve(controlNormalVectorsD.size());
    for (std::size_t i = 0; i < controlNormalOriginsD.size(); ++i) {
        const auto& p = controlNormalOriginsD[i];
        psControlNormalOriginsP.emplace_back(static_cast<float>(p[0]),
                                             static_cast<float>(p[1]),
                                             static_cast<float>(p[2]));
    }
    for (std::size_t i = 0; i < controlNormalVectorsD.size(); ++i) {
        const auto& v = controlNormalVectorsD[i];
        psControlNormalVectors.emplace_back(static_cast<float>(v[0]),
                                            static_cast<float>(v[1]),
                                            static_cast<float>(v[2]));
    }

    const auto controlsD = Utils::buildControlCloud(PF->getNeutralCurvenet());
    psControlsP.clear();
    psControlsP.reserve(controlsD.size());
    for (const auto& p : controlsD) {
        psControlsP.emplace_back(static_cast<float>(p[0]),
                                 static_cast<float>(p[1]),
                                 static_cast<float>(p[2]));
    }

    psCurvenet = polyscope::registerCurveNetwork("Curvenet", psCurvenetP, psCurvenetE);
    psCurvenet->setColor({0.95f, 0.55f, 0.15f});
    psCurvenet->setRadius(0.0018f, true);

    psDCurvenet = polyscope::registerCurveNetwork("DCurveNet Segments", psDCurvenetP, psDCurvenetE);
    psDCurvenet->setColor({0.15f, 0.80f, 0.95f});
    psDCurvenet->setRadius(0.0045f, true);

    std::vector<glm::vec3> plusT;
    std::vector<glm::vec3> plusB;
    std::vector<glm::vec3> plusN;
    std::vector<glm::vec3> minusT;
    std::vector<glm::vec3> minusB;
    std::vector<glm::vec3> minusN;
    const auto& dSegs = PF->getNeutralDCurvenet().segments();
    plusT.reserve(dSegs.size());
    plusB.reserve(dSegs.size());
    plusN.reserve(dSegs.size());
    minusT.reserve(dSegs.size());
    minusB.reserve(dSegs.size());
    minusN.reserve(dSegs.size());

    for (const auto& seg : dSegs) {
        if (seg.plusSide().valid) {
            plusT.push_back(toGlm(seg.plusSide().BS.col(0)));
            plusB.push_back(toGlm(seg.plusSide().BS.col(1)));
            plusN.push_back(toGlm(seg.plusSide().BS.col(2)));
        } else {
            plusT.emplace_back(0.f, 0.f, 0.f);
            plusB.emplace_back(0.f, 0.f, 0.f);
            plusN.emplace_back(0.f, 0.f, 0.f);
        }

        if (seg.minusSide().valid) {
            minusT.push_back(toGlm(seg.minusSide().BS.col(0)));
            minusB.push_back(toGlm(seg.minusSide().BS.col(1)));
            minusN.push_back(toGlm(seg.minusSide().BS.col(2)));
        } else {
            minusT.emplace_back(0.f, 0.f, 0.f);
            minusB.emplace_back(0.f, 0.f, 0.f);
            minusN.emplace_back(0.f, 0.f, 0.f);
        }
    }

    constexpr double kFrameDisplayScale = 1.0; // true geometric length

    auto* qPlusT = psDCurvenet->addEdgeVectorQuantity("Frame+ t (true)", plusT, polyscope::VectorType::STANDARD);
    qPlusT->setVectorColor({0.98f, 0.65f, 0.10f});
    qPlusT->setVectorLengthRange(1.0);
    qPlusT->setVectorLengthScale(kFrameDisplayScale, false);
    qPlusT->setEnabled(true);

    auto* qPlusB = psDCurvenet->addEdgeVectorQuantity("Frame+ b (true)", plusB, polyscope::VectorType::STANDARD);
    qPlusB->setVectorColor({0.95f, 0.20f, 0.45f});
    qPlusB->setVectorLengthRange(1.0);
    qPlusB->setVectorLengthScale(kFrameDisplayScale, false);
    qPlusB->setEnabled(true);

    auto* qPlusN = psDCurvenet->addEdgeVectorQuantity("Frame+ n (true)", plusN, polyscope::VectorType::STANDARD);
    qPlusN->setVectorColor({0.15f, 0.95f, 0.30f});
    qPlusN->setVectorLengthRange(1.0);
    qPlusN->setVectorLengthScale(kFrameDisplayScale, false);
    qPlusN->setEnabled(true);

    auto* qMinusT = psDCurvenet->addEdgeVectorQuantity("Frame- t (true)", minusT, polyscope::VectorType::STANDARD);
    qMinusT->setVectorColor({0.85f, 0.45f, 0.10f});
    qMinusT->setVectorLengthRange(1.0);
    qMinusT->setVectorLengthScale(kFrameDisplayScale, false);
    qMinusT->setEnabled(true);

    auto* qMinusB = psDCurvenet->addEdgeVectorQuantity("Frame- b (true)", minusB, polyscope::VectorType::STANDARD);
    qMinusB->setVectorColor({0.55f, 0.40f, 1.00f});
    qMinusB->setVectorLengthRange(1.0);
    qMinusB->setVectorLengthScale(kFrameDisplayScale, false);
    qMinusB->setEnabled(true);

    auto* qMinusN = psDCurvenet->addEdgeVectorQuantity("Frame- n (true)", minusN, polyscope::VectorType::STANDARD);
    qMinusN->setVectorColor({0.10f, 0.80f, 1.00f});
    qMinusN->setVectorLengthRange(1.0);
    qMinusN->setVectorLengthScale(kFrameDisplayScale, false);
    qMinusN->setEnabled(true);

    psDCurvenetInterior = polyscope::registerPointCloud("DCurveNet Vertices (Interior)", psDCurvenetInteriorP);
    psDCurvenetInterior->setPointColor({0.65f, 0.90f, 1.00f});
    psDCurvenetInterior->setPointRadius(0.0035, true);

    psDCurvenetControls = polyscope::registerPointCloud("DCurveNet Vertices (Controls)", psDCurvenetControlP);
    psDCurvenetControls->setPointColor({1.00f, 0.40f, 0.10f});
    psDCurvenetControls->setPointRadius(0.0075, true);

    psProjectedVertexSamples = polyscope::registerPointCloud("Projected Samples (Vertex)", psProjectedVertexSamplesP);
    psProjectedVertexSamples->setPointColor({1.00f, 0.25f, 0.20f});
    psProjectedVertexSamples->setPointRadius(0.0090, true);

    psProjectedEdgeSamples = polyscope::registerPointCloud("Projected Samples (Edge)", psProjectedEdgeSamplesP);
    psProjectedEdgeSamples->setPointColor({1.00f, 0.85f, 0.20f});
    psProjectedEdgeSamples->setPointRadius(0.0080, true);

    psProjectedFaceSamples = polyscope::registerPointCloud("Projected Samples (Face)", psProjectedFaceSamplesP);
    psProjectedFaceSamples->setPointColor({0.15f, 0.95f, 0.45f});
    psProjectedFaceSamples->setPointRadius(0.0070, true);

    psControls = polyscope::registerPointCloud("Controls", psControlsP);
    psControls->setPointRadius(0.0045, true);
    psControls->setPointColor({0.95f, 0.75f, 0.25f});

    psControlNormals = polyscope::registerPointCloud("Control Corner Normal Origins", psControlNormalOriginsP);
    psControlNormals->setPointColor({0.10f, 0.95f, 0.25f});
    psControlNormals->setPointRadius(0.005, true);

    if (!psControlNormalVectors.empty()) {
        constexpr double kControlNormalDisplayScale = 0.08; // visualization-only multiplier
        auto* nQ = psControlNormals->addVectorQuantity("Control Corner Normals",
                                                       psControlNormalVectors,
                                                       polyscope::VectorType::STANDARD);
        nQ->setEnabled(true);
        nQ->setVectorColor({0.10f, 0.95f, 0.25f});
        nQ->setVectorLengthRange(1.0);
        nQ->setVectorLengthScale(kControlNormalDisplayScale, false);
    }
}

// Performs call to pre-computation of cut-mesh and operators
int computePrecomp() {
    if (!PF) {
        std::cout << "Pre-computation skipped: profileformer not initialized." << std::endl;
        return 1;
    }
    PF->precomputation();
    refreshCurvenetVisuals();
    return 0;
}

// Performs call to surface deformation and updates PS mesh
int computeDeformation() {
    // TODO: Keep legacy entrypoint; deformation pipeline is not fully wired yet.
    return 0;
}

// Legacy aliases to preserve callback naming used by existing scaffold.
int performPrecomp() { return computePrecomp(); }
void endVertexEdit() { removeGizmo(); }
void resetVertexPositions() { resetMeshVertexPositions(); }


// Creates gizmo at vertex
void addGizmoAtVertex(int vert_idx) {
    const Eigen::Vector3d& p = psV[static_cast<std::size_t>(vert_idx)];
    glm::vec3 startpos(static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()));

    activeGizmo = true;
    if (!vertexGizmo) {
    vertexGizmo = polyscope::addTransformationGizmo("vertex_editor");
    vertexGizmo->setAllowTranslation(true);
    vertexGizmo->setAllowRotation(false);
    vertexGizmo->setAllowScaling(false);
    vertexGizmo->setInteractInLocalSpace(false);
    //vertexGizmo->setGizmoSize(0.5f);
    }

    // place gizmo at vertex position
    vertexGizmo->setPosition(startpos);
    return;
}

// Removes gizmo 
void removeGizmo() {
    if (!activeGizmo) return;

    if (vertexGizmo) {
        vertexGizmo->remove();
        vertexGizmo = nullptr;
    }
    activeGizmo = false;
    return;
}

// Modify a vertex position in polyscope
void modifyVertexPositions(int vert_idx, glm::vec3 new_pos) {
    psV[static_cast<std::size_t>(vert_idx)] =
        Eigen::Vector3d(static_cast<double>(new_pos.x),
                        static_cast<double>(new_pos.y),
                        static_cast<double>(new_pos.z));
    // Need a conversion function to change internal copy
    // V_current[vert_idx] = IO::glmToEigen(new_pos);
    return;
}

// Reset all vertex positions in Polyscope, as well as the constraints
/*
TODO: Some way to reset curvenet + controls
*/
void resetMeshVertexPositions() {
    // Reset Mesh
    Utils::copyPositions(V, psV);
    Utils::copyConnectivity(T, psT);
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", psV, psT);
    return;
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    ImGuiIO& io = ImGui::GetIO();
    bool mouseClicked = ImGui::IsMouseClicked(0);
    glm::vec2 screen{io.MousePos.x, io.MousePos.y};

    polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen);

    if (PF) {
        int samplesPerMeanEdge = PF->getDCurveSamplesPerMeanEdge();
        int uniformRefineSamples = PF->getDCurveUniformRefineSamples();
        bool discretizationChanged = false;

        ImGui::SeparatorText("Discrete Sampling");

        if (ImGui::SliderInt("Samples / Mean Edge", &samplesPerMeanEdge, 1, 32)) {
            discretizationChanged = true;
        }
        if (ImGui::SliderInt("Uniform Refine Samples", &uniformRefineSamples, 8, 512)) {
            discretizationChanged = true;
        }

        if (discretizationChanged) {
            PF->setDiscretizationParameters(samplesPerMeanEdge, uniformRefineSamples);
            computePrecomp();
        }

        const auto& dcn = PF->getNeutralDCurvenet();
        ImGui::Text("dVerts: %zu", dcn.verts().size());
        ImGui::Text("dSegments: %zu", dcn.segments().size());
    }

    if (ImGui::Button("Perform Pre-Computation")) {
        performPrecomp();
    }
    // Deformation stuff
    if (ImGui::Button("Compute Deformation")) {
        computeDeformation();
    }

    if (ImGui::Button(gizmoMode ? "Stop Moving Vertex" : "Move Vertex")) {
        gizmoMode = !gizmoMode;
        selectedVertex = -1;
        endVertexEdit();
    }

    ImGui::SameLine();
    if (ImGui::Button("Reset Vertices")) {
        std::cout << "Resetting vertex locations. Note that constraints will NOT be reset." << std::endl;
        resetVertexPositions();
        gizmoMode = false;
    }

    // User can apply controls --> These need to be in pairs and user must be able
    // to select previously selected controls
    if (ImGui::Button(controlMode ? "Stop Applying Controls" : "Apply Controls")) {
        controlMode = !controlMode;
    }

    // Constraint mode activation
    if (controlMode && mouseClicked && pick.isHit && pick.structure == psMesh) {
        polyscope::SurfaceMeshPickResult meshPick = psMesh->interpretPickResult(pick);

        // TODO: We should be able to click the mesh arbitrarily, not just at vertices
        if (meshPick.elementType == polyscope::MeshElement::VERTEX) {
            selectedVertex = static_cast<int>(meshPick.index);
        } else {
            std::cout << "Did not click on mesh." << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Too few arguments. Usage: ./profile_former <OBJ file path>" << std::endl;
        return 1;
    }
    InputPath = argv[1];
    CurvesPath = (argc > 2) ? argv[2] : "../curvenet/data/sphere-curves.json";

    // Initialize polyscope
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None; // Disable ground plane

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    polyscope::init();

    // Set camera view
    // polyscope::view::setUpDir(polyscope::UpDir::ZUp);      // Z up
    // polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront); // -Y forward

    // Set projection to orthographic
    // polyscope::view::setProjectionMode(polyscope::ProjectionMode::Orthographic);

    // Load our mesh object
    std::cout << "\nLoading surface mesh file" << std::endl;
    // Legacy line preserved (disabled): igl::readOBJ(InputPath, V, T);
    Utils::loadObjMesh(InputPath, V, T);

    std::cout << "Loading curve JSON file" << std::endl;
    JSONUtils::CurvenetInput cnInput = JSONUtils::loadBezierCurvenetInput(CurvesPath);

    // Initialize profileformer + curvenet with existing project types.
    PF = std::make_unique<ProfileFormer::profileformer>(V, T, cnInput.controlP, cnInput.surfaceN, cnInput.curveC);
    computePrecomp();

    // Register tet mesh with PS
    std::cout << "Registering Surface Mesh to Polyscope" << std::endl;
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", V, T);

    // Curvenet + controls were populated by computePrecomp() / refreshCurvenetVisuals()
    // // Empty curvenet
    // psCurvenet = polyscope::registerCurveNetwork("Discrete Curvenet", psCurvenetP, psCurvenetE);
    // //psCurvenet->setEnabled(true);

    // // Empty controls
    // psControls = polyscope::registerPointCloud("Controls", psCurvenetP);
    // //psControls->setEnabled(true);

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
