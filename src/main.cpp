#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <string>
#include <cmath>
#include <tuple>
#include <array>
#include <vector>
#include <stdexcept>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

//#include "args/args.hxx"
#include "imgui.h"

// My files
#include "profilemover/profilemover.hpp"
#include "psCurvenet/pscurvenet.hpp"
#include "utils/utils.hpp"
#include "IO/io.hpp"

// Main file for visualization with Polyscope

// VARIABLES FOR POLYSCOPE OBJECTS

// SURFACE MESH (M)
// NOTE: this needs to be a generic surface mesh and NOT a triangle mesh
Eigen::MatrixXd psMesh_V; // Vertex list
std::vector<std::vector<int>> psMesh_F; // Face list: Note the inner list has arbitrary size for non-triangle faces
polyscope::SurfaceMesh* psMesh = nullptr;

// CURVENET (CN)
Eigen::MatrixXd psCN_P; // Point list
std::vector<std::array<int, 2>> psCN_E; // Edge List
polyscope::CurveNetwork* psEditableCN = nullptr;

// CONTROLS (PC)
// Not sure what we need visualization-wise
Eigen::MatrixXd psControls_P; // Controls
polyscope::PointCloud* psControlsPC = nullptr;    // point cloud for controls
// TANGENTS (PC)
Eigen::MatrixXd psTangents_P; // Tangents
polyscope::PointCloud* psTangentsPC = nullptr;    // point cloud for tangents

// TANGENTS (CN)
Eigen::MatrixXd psTangentsVec; // Aggregate list of controls and tangents
std::vector<std::array<int, 2>> psTangents_E; // Edge List between controls and tangents
polyscope::CurveNetwork* psTangentsCN = nullptr;  // Connects controls to their tangents


// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;
std::string OutputPath;

// UI HELPERS
bool precompDone = false;   // Once precomputation is done, we can no longer edit the curves

bool createCtrlMode = false;  // Allows users to place control points
bool createSplineMode = false;  // Allow users to initialize new splines

bool editCtrlMode = false;   // Allows users to modify controls
bool editTanMode = false;   // Allows users to modify tangents

bool delCtrlMode = false;   // Allows users to remove control points
bool delSplineMode = false;  // Allows users to remove splines

// Spline creation/removal helpers
int selectedIdx = -1;    // Index of selected vertex on mesh
std::pair<int, int> selectedPair = {-1, -1};

// Editing helpers
bool tanConstraint = true;  // Constrain tangent movement to tangent plane only

// Gizmo helpers
bool activeGizmo = false; // This tells us if there is an active gizmo
Eigen::Vector3d gizmoPos;
static polyscope::TransformationGizmo* vertexGizmo = nullptr;

// Pre-computation
int samplingParam = 5;

// Discrete curvenet for modeling
std::unique_ptr<psCurvenet::pscurvenet> psCN = nullptr; // Curvenet that polyscope will use for updates

// Profile Mover
//ProfileMover::profilemover PM;


// ----------------- FUNCTIONS BEGIN HERE -------------------------

// Performs call to pre-computation of cut-mesh and operators
int performPrecomp() {
    // Make a hard copy of the neutral curvenet and pass it as the profile mover's copy
    // That way we maintain one copy to that polyscope can edit
    // TODO: Initialize global profile mover object

    // TODO: Time the pre-computation to see how time-intensive it is.
    // TODO: Catch and handle errors here
    return 1;
}

// Performs call to surface deformation and updates PS mesh
int computeDeformation() {
    return 1;
}

// Saves current curvenet to some file format
int saveCurvenet() {
    std::cout << "Curvenet saved to file." << std::endl;
    return 1;
}

// Creates gizmo at vertex
void addGizmoAtLocation(Eigen::Vector3d& startpos) {
    activeGizmo = true;
    if (!vertexGizmo) {
        vertexGizmo = polyscope::addTransformationGizmo("vertex_editor");
        vertexGizmo->setAllowTranslation(true);
        vertexGizmo->setAllowRotation(false);
        vertexGizmo->setAllowScaling(false);
        vertexGizmo->setInteractInLocalSpace(false);
        //vertexGizmo->setGizmoSize(0.5f);
    }

    // Place gizmo at vertex position
    vertexGizmo->setPosition(Utils::eigenToGLM(startpos));
    return;
}
// Creates a new gizmo rotated to match an axis
void addGizmoAtLocation(Eigen::Vector3d& startpos, Eigen::Vector3d& ax1, Eigen::Vector3d& ax2, Eigen::Vector3d& ax3) {
    activeGizmo = true;
    if (!vertexGizmo) {
        glm::mat4 T(1.0f);
        T[0] = glm::vec4(Utils::eigenToGLM(ax1), 0.0f);
        T[1] = glm::vec4(Utils::eigenToGLM(ax2), 0.0f);
        T[2] = glm::vec4(Utils::eigenToGLM(ax3), 0.0f);
        T[3] = glm::vec4(Utils::eigenToGLM(startpos), 1.0f);
        vertexGizmo = polyscope::addTransformationGizmo("editor");
        vertexGizmo->setTransform(T);
        vertexGizmo->setAllowTranslation(true);
        vertexGizmo->setAllowRotation(true);
        vertexGizmo->setAllowScaling(false);
        vertexGizmo->setInteractInLocalSpace(true); // TODO: Check this
        //vertexGizmo->setGizmoSize(0.5f);
    }

    // place gizmo at vertex position
    vertexGizmo->setPosition(Utils::eigenToGLM(startpos));
    gizmoPos = startpos;
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

// Remove all polyscope object
void removeAllCurvenetPS() {
    if (psEditableCN) {
        psEditableCN->remove();
        psEditableCN = nullptr;
    }

    if (psTangentsCN) {
        psTangentsCN->remove();
        psTangentsCN = nullptr;
    }

    if (psControlsPC) {
        psControlsPC->remove();
        psControlsPC = nullptr;
    }

    if (psTangentsPC) {
        psTangentsPC->remove();
        psTangentsPC = nullptr;
    }
    return;
}

void updateCurvenet(bool conn = false) {
    // Reset curvenet
    psCN->cnAsCurveNetwork(psCN_P, psCN_E);
    psCN->tansAsCurveNetwork(psTangentsVec, psTangents_E);
    psCN->cPosAsMatrix(psControls_P);
    psCN->tPosAsMatrix(psTangents_P);

    if (conn) {
        removeAllCurvenetPS();
        if (psCN_E.size() > 0) {
            psEditableCN = polyscope::registerCurveNetwork("Curvenet", psCN_P, psCN_E);
            psEditableCN->setColor({0.0f, 0.0f, 0.0f});
            psEditableCN->setMaterial("flat");
            psEditableCN->setTransparency(0.65);
            psEditableCN->setRadius(0.003);
            psEditableCN->setEnabled(true);
        }

        if (psTangents_E.size() > 0) {
            psTangentsCN = polyscope::registerCurveNetwork("Handles", psTangentsVec, psTangents_E);
            psTangentsCN->setColor({0.5f, 0.55f, 0.15f});
            psTangentsCN->setMaterial("flat");
            psTangentsCN->setTransparency(0.8);
            psTangentsCN->setRadius(0.006);
            psTangentsCN->setEnabled(true);
        }

        if (psControls_P.rows() > 0) {
            psControlsPC = polyscope::registerPointCloud("Controls", psControls_P);
            psControlsPC->setPointColor({0.9f, 0.2f, 0.1f});
            psControlsPC->setMaterial("flat");
            psControlsPC->setPointRadius(0.02);
            psControlsPC->setEnabled(true);
        }

        if (psTangents_P.rows() > 0) {
            psTangentsPC = polyscope::registerPointCloud("Tangents", psTangents_P);
            psTangentsPC->setPointColor({0.1f, 0.9f, 0.2f});
            psTangentsPC->setMaterial("flat");
            psTangentsPC->setPointRadius(0.012);
            psTangentsPC->setEnabled(true);
        }
    } else {
        if (psEditableCN && psCN_E.size() > 0) {
            psEditableCN->updateNodePositions(psCN_P);
        }
        if (psTangentsCN && psTangents_E.size() > 0) {
            psTangentsCN->updateNodePositions(psTangentsVec);
        }
        if (psControlsPC && psControls_P.rows() > 0) {
            psControlsPC->updatePointPositions(psControls_P);
        }
        if (psTangentsPC && psTangents_P.rows() > 0) {
            psTangentsPC->updatePointPositions(psTangents_P);
        }
    }
    return;
}

// Reset all control positions in Polyscope
void resetCurvenet() {
    // If not precomputed
    if (!precompDone) {
        std::cout << "Cannot reset. Precomputation not performed." << std::endl;
        return;
    }
    // TODO: Reset curvenet
    return;
}

// clear all modes and their variables except for the specified mode
int clearModes() {
    createCtrlMode = false;
    createSplineMode = false;
    editCtrlMode = false;
    editTanMode = false;
    delCtrlMode = false;
    delSplineMode = false;
    
    removeGizmo();

    selectedIdx = -1;
    gizmoPos = Eigen::Vector3d::Zero();

    selectedPair = {-1, -1};
    return 1;
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    ImGuiIO& io = ImGui::GetIO();
    bool mouseDown = ImGui::IsMouseDown(0);
    bool mouseClicked = ImGui::IsMouseClicked(0);
    bool mouseReleased = ImGui::IsMouseReleased(0);
    glm::vec2 screen{io.MousePos.x, io.MousePos.y};

    polyscope::PickResult pick = polyscope::pickAtScreenCoords(screen);

    // Pre-compute cut-mesh and operators
    if (ImGui::Button("Perform Pre-Computation")) {
        clearModes();
        if (precompDone) {
            std::cout << "Pre-computation on neutral pose already performed." << std::endl;
        } else if (psControls_P.rows() <= 1) {
            std::cout << "No splines specified. Add one or more spline before pre-computing." << std::endl;
        } else {
            // Compress the curvenet
            psCN->cleanupControls();
            updateCurvenet();
            // Precompute
            std::cout << "Performing pre-computation on neutral pose." << std::endl;
            performPrecomp();
            precompDone = true;
        }
    }

    // Deformation stuff
    // TODO: In the future, make this mode automatic after running precomp
    if (ImGui::Button("Deformation Mode")) {
        clearModes();
        if (!precompDone) {
            std::cout << "Perform pre-computation before applying deformation." << std::endl;
        } else {
            std::cout << "Computing Deformation." << std::endl;
            computeDeformation();
        }
    }

    // User parameter for sampling the spline
    ImGui::SliderInt("Sampling Param", &samplingParam, 2, 8);

    // Save the current curvenet state
    if (ImGui::Button("Save Curvenet")) {
        clearModes();
        saveCurvenet();
    }


    // CONTROL/SPLINE CREATION
    // Create controls
    if (ImGui::Button(createCtrlMode ? "Stop Creating Controls" : "Create Controls")) {
        bool tempMode = createCtrlMode;
        clearModes();
        createCtrlMode = !tempMode;
    }
    ImGui::SameLine();
    // Move a curvenet vertex
    if (ImGui::Button(createSplineMode ? "Stop Creating Splines" : "Create Splines")) {
        bool tempMode = createSplineMode;
        clearModes();
        createSplineMode = !tempMode;
    }


    // CONTROL/TANGENT EDITING
    if (ImGui::Button(editCtrlMode ? "Stop Moving Control" : "Move Control")) {
        bool tempMode = editCtrlMode;
        clearModes();
        if (psControls_P.rows() == 0) {
            std::cout << "Cannot edit. No existing controls." << std::endl;
            editCtrlMode = false;
        } else if (tempMode == false) {
            editCtrlMode = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button(editTanMode ? "Stop Moving Splines" : "Move Splines")) {
        bool tempMode = editTanMode;
        clearModes();
        if (psTangents_P.rows() == 0) {
            std::cout << "Cannot edit. No existing splines." << std::endl;
            editTanMode = !tempMode;
        } else if (tempMode == false) {
            editTanMode = true;
        }
    }
    ImGui::SameLine();
    ImGui::Checkbox("Proj. Tans", &tanConstraint); // TODO: do not allow degenerate vectors --> Constrain tan vertex AND gizmo

    // CONTROL/SPLINE DELETION
    if (ImGui::Button(delCtrlMode ? "Stop Removing Controls" : "Remove Controls")) {
        bool tempMode = delCtrlMode;
        clearModes();
        if (psControls_P.rows() == 0) {
            std::cout << "Cannot delete. No existing controls." << std::endl;
            delCtrlMode = false;
        } else if (tempMode == false) {
            delCtrlMode = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button(delSplineMode ? "Stop Removing Splines" : "Remove Splines")) {
        bool tempMode = delSplineMode;
        clearModes();
        if (psTangents_P.rows() == 0) {
            std::cout << "Cannot delete. No existing splines." << std::endl;
            delSplineMode = !tempMode;
        } else if (tempMode == false) {
            delSplineMode = true;
        }
    }
    
    // RESETs
    if (ImGui::Button("Clear Gizmo")) {
        std::cout << "Removing current gizmo." << std::endl;
        removeGizmo();
        selectedIdx = -1;
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear Curvenet")) {
        std::cout << "Clearing entire curvenet." << std::endl;
        psCN->resetCurvenet();
        updateCurvenet(true);
    }
    // May need to store a copy of the rest curvenet
    if (ImGui::Button("Reset Curvenet")) {
        std::cout << "Resetting Splines." << std::endl;
        resetCurvenet();
        clearModes();
    }

    // CREATE MODE CLICKS
    // Clicked on mesh during create mode
    if (createCtrlMode && mouseClicked && pick.isHit) {
        // If we click a current control, intializes the control and tangents
        if (pick.structure == psMesh) {
            polyscope::SurfaceMeshPickResult meshPick = psMesh->interpretPickResult(pick);
            Eigen::Vector3d pos = Utils::glmToEigen(pick.position);

            // TODO: Switch to a different version for non-triangle meshes
            Eigen::Vector3d normal;
            int valid = Utils::closestPointNormalOnMesh(pos, psMesh_V, psMesh_F, normal);
            if (valid == 1) {
                psCN->addControl(pos, normal);
                updateCurvenet(true);
                std::cout << "New Vert created at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
            } else {
                std::cout << "No valid point picked." << std::endl;
            }
        }
    }
    // Create Spline by picking two controls
    if (createSplineMode && mouseClicked && pick.isHit) {
        if (pick.structure == psControlsPC) {
            polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);

            if (selectedPair.first == -1) {
                selectedPair.first = selectedIdx;
                Eigen::Vector3d pos = psControls_P.row(selectedIdx).transpose();
                std::cout << "First Spline Vert: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
            } else {
                selectedPair.second = selectedIdx;
                Eigen::Vector3d pos = psControls_P.row(selectedIdx).transpose();
                std::cout << "Second Spline Vert: (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
                // Compute initial tangent directions
                psCN->addSpline(selectedPair.first, selectedPair.second);
                // Reset pair
                selectedPair = {-1, -1};
                std::cout << "New Spline Created.\n" << std::endl;
                updateCurvenet(true);
            }
        }
    }

    // REMOVE MODE CLICKS
    // Clicked on control to remove
    if (delCtrlMode && mouseClicked && pick.isHit && pick.structure == psControlsPC) {
        polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

        psCN->removeControl(static_cast<int>(pcPick.index));
        std::cout << "Control removed." << std::endl;
        updateCurvenet(true);
    }
    // Clicked on a tangent whose spline we should remove
    if (delSplineMode && mouseClicked && pick.isHit && pick.structure == psTangentsPC) {
        polyscope::PointCloudPickResult pcPick = psTangentsPC->interpretPickResult(pick);

        psCN->removeSplineByTangent(static_cast<int>(pcPick.index));
        std::cout << "Spline removed." << std::endl;
        updateCurvenet(true);
    }

    // EDIT MODE CLICKS
    // Select control to edit
    if (editCtrlMode && mouseClicked) {
        if (pick.isHit && pick.structure == psControlsPC) {
            polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);
            // Build a basis
            Eigen::Vector3d selectedPos = psControls_P.row(selectedIdx).transpose();
            Eigen::Vector3d selectedN = psCN->getNormal(selectedIdx);
            Eigen::Vector3d t0, t1;
            Utils::buildPlaneBasis(selectedN, t0, t1);
            std::cout << "Editing Vert at (" << selectedPos[0] << ", " << selectedPos[1] << ", " << selectedPos[2] << ")" << std::endl;
            // add Gizmo at position
            addGizmoAtLocation(selectedPos, selectedN, t0, t1);
        } else if (editCtrlMode && mouseClicked && !activeGizmo && (!pick.isHit || (pick.isHit && pick.structure != psControlsPC))) {  // or clear
            clearModes();
            editCtrlMode = true;
        }
    }
    // Select tangent to edit
    if (editTanMode && mouseClicked) {
        if (pick.isHit && pick.structure == psTangentsPC) {
            // Index into tangent list. Get associated spline by integer dividing by 2.
            // Then use the spline index to edit the tangent's position directly
            polyscope::PointCloudPickResult pcPick = psTangentsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);
            Eigen::Vector3d selectedPos = psTangents_P.row(selectedIdx).transpose();
            std::cout << "Editing Spline Handle at (" << selectedPos[0] << ", " << selectedPos[1] << ", " << selectedPos[2] << ")" << std::endl;
            // add Gizmo at position
            addGizmoAtLocation(selectedPos);
        } else if (editTanMode && mouseClicked && !activeGizmo && (!pick.isHit || (pick.isHit && pick.structure != psTangentsPC))) {  // or clear
            clearModes();
            editTanMode = true;
        }
    }

    // Update control position per-frame
    if (editCtrlMode && activeGizmo && selectedIdx >= 0) {
        // Get the gizmo's location at this frame
        Eigen::Vector3d gizmoPosF = Utils::glmToEigen(vertexGizmo->getPosition());
        glm::mat4 T = vertexGizmo->getTransform();
        Eigen::Vector3d gizmoNormal = Utils::glmToEigen(glm::normalize(glm::vec3(T[0])));

        psCN->updateControlPos(selectedIdx, gizmoPosF, tanConstraint);
        psCN->updateControlNormal(selectedIdx, gizmoNormal, true, tanConstraint);
        // Update the associate tangent point in the curvenet
        vertexGizmo->setPosition(Utils::eigenToGLM(gizmoPosF));

        updateCurvenet();
    }
    // Update tangent position per-frame
    if (editTanMode && activeGizmo && selectedIdx >= 0) {
        // Get the gizmo's location at this frame
        Eigen::Vector3d gizmoPosF = Utils::glmToEigen(vertexGizmo->getPosition());
        // If we are too close to either endpoint, do not update
        bool updated = psCN->updateTangentPos(selectedIdx, gizmoPosF, tanConstraint);
        updateCurvenet();
        if (updated) {
            Eigen::Vector3d tangentPos = psTangents_P.row(selectedIdx).transpose();
            vertexGizmo->setPosition(Utils::eigenToGLM(tangentPos));
        }
    }

    return;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Too few arguments. Usage: ./profile_mover <input OBJ file path>" << std::endl;
        return 1;
    }
    InputPath = argv[1];
    // OutputPath = argv[2];

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
    IO::readOBJ(InputPath, psMesh_V, psMesh_F);

    // Register mesh with PS
    std::cout << "Registering Surface Mesh to Polyscope" << std::endl;
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", psMesh_V, psMesh_F);
    psMesh->setSurfaceColor({0.3f, 0.2f, 1.0f});

    // Create polyscope's curvenet
    psCN = std::make_unique<psCurvenet::pscurvenet>();

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}