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
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

//#include "args/args.hxx"
#include "imgui.h"

// My files
#include "profilemover/profilemover.hpp"
#include "curvenet/components/spline.hpp"   // So we can segment splines on the fly for visualization
#include "mesh/mesh.hpp"
#include "utils/utils.hpp"
#include "IO/io.hpp"

// Main file for visualization with Polyscope

/*
NOTES:
- We store two copies of each object. One is for Polyscope to mess with, the other is the "neutral"/original copy
- TODO: When a control moves, the associated tangents should all move by the same amount
        Write a function which converts a curvenet to a polyscope curve network/ point cloud object
*/

// VARIABLES FOR POLYSCOPE OBJECTS
// Handles for the surface mesh
Eigen::MatrixXd psV; // Vertex list
std::vector<std::vector<int>> psF; // Face list: Note the inner list has arbitrary size for non-triangle faces
polyscope::SurfaceMesh* psMesh;

// Handles for discrete curvenet
Eigen::MatrixXd psCurvenetP; // Point list
std::vector<std::array<int, 2>> psCurvenetE; // Edge List
polyscope::CurveNetwork* psCurvenet;

// Handles for controls (?)
// Not sure what we need visualization-wise
Eigen::MatrixXd psControlsP; // Controls
Eigen::MatrixXd psTangentsP; // Tangents
polyscope::PointCloud* psControlsPC;    // point cloud for controls
polyscope::PointCloud* psTangentsPC;    // point cloud for tangents

Eigen::MatrixXd psControlsAndTangents; // Aggregate list of controls and tangents
std::vector<std::array<int, 2>> psControlsAndTangentsE; // Edge List between controls and tangents
std::map<int, int> tanListToCTList; // Maps points in the tangent list to points in the psControlsAndTangents list
polyscope::CurveNetwork* psControlsCN;  // Connects controls to their tangents


// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;
std::string OutputPath;

// UI HELPERS
bool precompDone = false;   // Once precomputation is done, we can no longer edit the curves
bool gizmoMode = false; // This allows the user to create a gizmo
bool createMode = false;  // Allows users to place control points
bool tanMode = false;   // Allows users to modify controls and tangents
bool removeMode = false;   // Allows users to remove control points
bool recomputeCurvenet = false; // Whether we should recompute the ps curve networks

// Spline creation/removal helpers
int selectedIdx = -1;    // Index of selected vertex on mesh
int num_selected = 0;   // Number valid entries selected so far
std::array<Eigen::Vector3d, 2> selectedPair;
std::array<Eigen::Vector3d, 2> selectedPairNormals;
std::array<int, 2> selectedPairIdx;  // Flag vertices that come from the existing point cloud

// Editing helpers
bool editObjectSelected = false;    // Whether we have selected a vertex to add a gizmo to
int controlIdx = -1;    // Index of control selected OR control associated with tangent
int tangentIdx = -1;    // Index of tangent selected to be edited
bool tanConstraint = true;  // Constrain tangent movement to tangent plane only

// Gizmo helpers
bool activeGizmo = false; // This tells us if there is an active gizmo
Eigen::Vector3d gizmoPos;
static polyscope::TransformationGizmo* vertexGizmo = nullptr;

// Pre-computation
int samplingParam = 5;

// INTERNAL OBJECT COPIES
// Internal mesh copy (Eigen)
// Neutral mesh state
std::vector<Eigen::Vector3d> V;
std::vector<std::vector<int>> F;
std::unique_ptr<Mesh::mesh> neutralMesh;

// Neutral discrete curvenet
std::unique_ptr<Curvenet::curvenet> editableCurvenet; // Curvenet that polyscope will use for updates
Curvenet::curvenet neutralCurvenet; // Copy of neutral curvenet to be initialized during pre-computation

// Profile Mover
//ProfileMover::profilemover PM;


// ----------------- FUNCTIONS BEGIN HERE -------------------------

// Performs call to pre-computation of cut-mesh and operators
int performPrecomp() {
    // Make a hard copy of the neutral curvenet and pass it as the profile mover's copy
    // That way we maintain one copy to that polyscope can edit
    // TODO: Initialize global profile mover object

    // TODO: Time the pre-computation to see how time-intensive it is.
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

    // place gizmo at vertex position
    vertexGizmo->setPosition(Utils::eigenToGLM(startpos));
    return;
}

void addGizmoAtLocation(Eigen::Vector3d& startpos, Eigen::Vector3d& ax1, Eigen::Vector3d& ax2, Eigen::Vector3d& ax3) {
    activeGizmo = true;
    if (!vertexGizmo) {
        glm::mat4 T(1.0f);
        T[0] = glm::vec4(Utils::eigenToGLM(ax1), 0.0f);
        T[1] = glm::vec4(Utils::eigenToGLM(ax2), 0.0f);
        T[2] = glm::vec4(Utils::eigenToGLM(ax3), 0.0f);
        T[3] = glm::vec4(Utils::eigenToGLM(startpos), 1.0f);
        vertexGizmo = polyscope::addTransformationGizmo("vertex_editor");
        vertexGizmo->setTransform(T);
        vertexGizmo->setAllowTranslation(true);
        vertexGizmo->setAllowRotation(false);
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

// Modify a vertex position in polyscope
void modifyVertexPositions(glm::vec3 new_pos) {
    // TODO: This needs to be tangent or controls, NOT mesh psV
    // psV[vert_idx] = new_pos;
    // Need a conversion function to change internal copy
    return;
}

// Reset all control positions in Polyscope
/*
TODO: Do NOT let user reset mesh. Can only reset curvenet. Notice how resetting curvenet should reset mesh.
*/
void resetVertexPositions() {
    // Reset
    return;
}

// clear all modes and their variables except for the specified mode
int clearModes(int exception = -1, bool clear_selection = true) {
    if (exception != 1) {
        createMode = false;
    }
    if (exception != 2) {
        removeMode = false;
    }
    if (exception != 3) {
        tanMode = false;
    }

    if (clear_selection) {
        removeGizmo();

        editObjectSelected = false;
        controlIdx = -1;
        tangentIdx = -1;
        gizmoPos.fill(0.0);

        selectedIdx = -1;
        num_selected = 0;
        selectedPair.fill(Eigen::Vector3d::Zero());
        selectedPairNormals.fill(Eigen::Vector3d::Zero());
        selectedPairIdx.fill(-1);
    }
    return exception;
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
        } else if (psCurvenetP.size() <= 1) {
            std::cout << "No splines specified. Add one or more spline before pre-computing." << std::endl;
        } else {
            std::cout << "Performing pre-computation on neutral pose." << std::endl;
            performPrecomp();
            precompDone = true;
        }
    }

    // Deformation stuff
    if (ImGui::Button("Compute Deformation")) {
        clearModes();
        if (!precompDone) {
            std::cout << "Perform pre-computation before applying deformation." << std::endl;
        } else {
            std::cout << "Computing Deformation." << std::endl;
            computeDeformation();
        }
    }

    // User parameter for sampling the spline
    // ImGui::SliderInt("Sampling Param", &samplingParam, 1, 10);

    // Save the current curvenet state
    if (ImGui::Button("Save Curvenet")) {
        clearModes();
        saveCurvenet();
    }

    // Move a curvenet vertex
    if (ImGui::Button(gizmoMode ? "Stop Moving Vertex" : "Move Vertex")) {
        clearModes();
        if (!precompDone) {
            std::cout << "No pre-computation performed. Apply pre-computation first." << std::endl;
        }
        gizmoMode = !gizmoMode;
    }

    ImGui::SameLine();
    // TODO: Do NOT let user modify underlying mesh. Can only modify the controls
    // May need to store a copy of the rest curvenet
    if (ImGui::Button("Reset Vertices")) {
        std::cout << "Resetting Controls." << std::endl;
        resetVertexPositions();
        clearModes();
        recomputeCurvenet = true;
    }

    // User can apply controls --> These need to be in pairs and user must be able
    // to select previously selected controls
    if (ImGui::Button(createMode ? "Stop Creating Splines" : "Create Spline")) {
        clearModes(1);
        createMode = !createMode;
    }
    
    // NOTE: May want to remove this. Indexing becomes annoying.
    ImGui::SameLine();
    if (ImGui::Button(removeMode ? "Stop Removing Splines" : "Remove Spline")) {
        // TODO: Implement spline removal
        clearModes(2);
        removeMode = !removeMode;
    }

    // TODO: Before pre-comp, only allow editing of tangents. 
    if (ImGui::Button(tanMode ? "Stop Editing" : "Edit Handles")) {
        clearModes(3);
        tanMode = !tanMode;
    }
    ImGui::SameLine();
    ImGui::Checkbox("Constrain to Tan. Plane", &tanConstraint); // TODO: do not allow degenerate vectors --> Constrain tan vertex AND gizmo


    // CREATE MODE CLICKS

    // Clicked on mesh during create mode
    // If use clicked away after selecting the first vertex, then clear it
    /*if (createMode && mouseClicked && (!pick.isHit || (pick.structure != psMesh && pick.structure != psControlsPC))) {
        if (num_selected == 1) {
            std::cout << "User clicked away and pair has been cleared. Please select your first control point." << std::endl;
            clearModes(1);
        }
    }*/

    // Clicked on mesh during create mode
    if (createMode && mouseClicked && pick.isHit) {
        if (pick.structure == psControlsPC) {
            polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

            selectedIdx = static_cast<int>(pcPick.index);
            Eigen::Vector3d pick_pos = psControlsP.row(selectedIdx);
            selectedPair[num_selected] = pick_pos;
            selectedPairIdx[num_selected] = selectedIdx;
            selectedPairNormals[num_selected] = (editableCurvenet->controlPoints)[selectedIdx].getNormal();
            num_selected++;
        } else if (pick.structure == psMesh) {
            polyscope::SurfaceMeshPickResult meshPick = psMesh->interpretPickResult(pick);
            // Check what mesh element type we hit
            if (meshPick.elementType == polyscope::MeshElement::VERTEX) {
                selectedIdx = static_cast<int>(meshPick.index);
                Eigen::Vector3d pick_pos = psV.row(selectedIdx);
                selectedPair[num_selected] = pick_pos;
                selectedPairIdx[num_selected] = -1;
                // Get vertex normal
                selectedPairNormals[num_selected] = neutralMesh->getVNormal(selectedIdx);
                num_selected++;
            } else if (meshPick.elementType == polyscope::MeshElement::FACE) {
                // Special exception for edges: need to compute nearest edge on the selected face
                // NOTE: Since we don't have a true edge-picker, check if face-pick is within some threshold
                Eigen::Vector3d pick_pos = Utils::glmToEigen(pick.position);
                Eigen::Vector3d nearestPoint;
                int tempIdx;
                double dist = neutralMesh->computeNearestFaceEdge(static_cast<int>(meshPick.index),
                                                                    pick_pos, tempIdx, nearestPoint);
                selectedPairIdx[num_selected] = -1;
                if (dist <= 1e-8 * neutralMesh->getMeanE()) {    // if sufficiently close to an edge
                    selectedPair[num_selected] = nearestPoint;
                    // Get edge normal
                    selectedPairNormals[num_selected] = neutralMesh->getENormal(tempIdx);
                } else {
                    selectedIdx = static_cast<int>(meshPick.index);
                    selectedPair[num_selected] = pick_pos;
                    // Get face normal
                    selectedPairNormals[num_selected] = neutralMesh->getFNormal(selectedIdx);
                }
                num_selected++;
            } else {
                std::cout << "Selected non- face or vertex part of mesh. No action performed." << std::endl;
            }
            /*
        // NOTE: Polyscope does not properly support edge-picking for non-triangular faces
        // so for safety I'm turning this off.
        else if (meshPick.elementType == polyscope::MeshElement::EDGE) {
            // Get edge normal
            selectedPairNormals.push_back(neutralMesh->getENormal(selectedIdx));
        }
        */
            if (num_selected == 1) {
                std::cout << "1st Point Selected: (" << selectedPair[0](0) << ", "
                                                   << selectedPair[0](1) << ", " 
                                                   << selectedPair[0](2) << ")" << std::endl;
            }
        } else {
            std::cout << "Did not click on the mesh or an existing control. Please try again." << std::endl;
        }
    }

    // If we have selected two verts to create, add them to spline
    if (createMode && (num_selected == 2)) {
        std::cout << "2nd Point Selected: (" << selectedPair[1](0) << ", "
                                            << selectedPair[1](1) << ", "
                                            << selectedPair[1](2) << "), " << std::endl;
        // Create tangents by finding vector and adding it to point
        Eigen::Vector3d toVec = selectedPair[1] - selectedPair[0];
        double scale = toVec.norm() * (1.2e-1) * (neutralMesh->getBBoxDiag());
        Eigen::Vector3d t0;
        Utils::projectVectorOntoTangentPlane(selectedPairNormals[0], toVec, t0, scale);
        t0 += selectedPair[0];
        Eigen::Vector3d t1;
        Utils::projectVectorOntoTangentPlane(selectedPairNormals[1], -1 * toVec, t1, scale);
        t1 += selectedPair[1];
        std::array<Eigen::Vector3d, 2> addedTangents = {t0, t1};
        std::cout << "Tangent points: (" << addedTangents[0](0) << ", "
                                            << addedTangents[0](1) << ", "
                                            << addedTangents[0](2) << "), ("
                                            << addedTangents[1](0) << ", "
                                            << addedTangents[1](1) << ", "
                                            << addedTangents[1](2) << ")" << std::endl;
        
        // Add spline to curvenet
        
        // Check if we have any existing control points
        if ((selectedPairIdx[0] != -1) && (selectedPairIdx[1] != -1)) { // Both exist
            editableCurvenet->addSpline(selectedPairIdx, addedTangents);
        } else if ((selectedPairIdx[0] == -1) && (selectedPairIdx[1] == -1)) { // Neither exists
            std::array<Eigen::Vector3d, 2> addedControls = {selectedPair[0], selectedPair[1]};
            std::array<Eigen::Vector3d, 2> addedNormals = {selectedPairNormals[0], selectedPairNormals[1]};
            editableCurvenet->addSpline(addedControls, addedNormals, addedTangents);
        } else if (selectedPairIdx[0] == -1) {  // Only the second one exists
            editableCurvenet->addSpline(selectedPair[0], selectedPairIdx[1], selectedPairNormals[0], addedTangents);
        } else if (selectedPairIdx[1] == -1) {  // Only the first one exists
            editableCurvenet->addSpline(selectedPairIdx[0], selectedPair[1], selectedPairNormals[1], addedTangents);
        }
        
        recomputeCurvenet = true;
        // clear relevant variables
        clearModes();
    }

    // REMOVE MODE CLICKS

    // Clicked on point during removeMode
    if (removeMode && mouseClicked && pick.isHit && pick.structure == psControlsPC) {
        polyscope::PointCloudPickResult pcPick = psControlsPC->interpretPickResult(pick);

        // selectedIdx = static_cast<int>(pcPick.index);

        // If we already selected one previously, then compute tangents
        
            // Check if such a spline exists between these two points

            // If so,
            // Remove the spline segments going between them + the tangents
            // Check if either has any splines anymore
            // If not, remove that control point and its tangents
        // Otherwise, mark for removal
    }

    // EDIT MODE CLICKS
    // TODO: Note that after pre-computation, we need to apply deformation
    // For pre-computation, we just need to change the tangent locations in the curvenet. How do we index...?
    if (tanMode && mouseClicked && pick.isHit && pick.structure == psTangentsPC) {
        editObjectSelected = true;
        // Index into tangent list. Get associated spline by integer dividing by 2.
        // Then use the spline index to edit the tangent's position directly
        polyscope::PointCloudPickResult pcPick = psTangentsPC->interpretPickResult(pick);

        tangentIdx = static_cast<int>(pcPick.index);
        Eigen::Vector3d pick_pos = psTangentsP.row(tangentIdx);
        int splineIdx = tangentIdx / 2;
        // Get parent control index
        if ((editableCurvenet->splines)[splineIdx].t0 == tangentIdx) {
            controlIdx = (editableCurvenet->splines)[splineIdx].c0;
        } else {
            controlIdx = (editableCurvenet->splines)[splineIdx].c1;
        }
        Eigen::Vector3d tanPos = psTangentsP.row(tangentIdx);
        // add Gizmo at position
        addGizmoAtLocation(tanPos);
    } else if (tanMode && mouseClicked && pick.isHit && pick.structure != psTangentsPC) {
        editObjectSelected = false;
        removeGizmo();
        clearModes(3);
    }
    // At each iter., update the position
    if (tanMode && activeGizmo) {
        // Get the gizmo's location at this frame
        Eigen::Vector3d gizmoPosF = Utils::glmToEigen(vertexGizmo->getPosition());
        Eigen::Vector3d controlPos = psControlsP.row(controlIdx);
        if (tanConstraint) {
            gizmoPosF = Utils::projectPointOntoPlane((editableCurvenet->controlPoints)[controlIdx].getNormal(), 
                                                    (editableCurvenet->controlPoints)[controlIdx].getPos(),
                                                    gizmoPosF);
        }
        // If we are too close to normal, copy in the old position.
        if ((gizmoPosF - controlPos).norm() <= 5e-2) {
            gizmoPosF = gizmoPos;
        } else {    // Otherwise, update our current gizmo position
            gizmoPos = gizmoPosF;
        }
        // Update the associate tangent point in the curvenet
        vertexGizmo->setPosition(Utils::eigenToGLM(gizmoPosF));
        (editableCurvenet->tangentPoints)[tangentIdx].setPos(gizmoPosF);

        recomputeCurvenet = true;
    }


    // Recompute the curvenet based on the updated information
    if (recomputeCurvenet) {
        // Recompute sampled spline and update polyscope
        // TODO: set vertex and edge sizes for various objects
        std::vector<Eigen::Vector3d> P;
        editableCurvenet->convertCurvnetToCN(P, psCurvenetE);
        IO::convertVertsToMatrix(P, psCurvenetP);
        psCurvenet = polyscope::registerCurveNetwork("Curvenet", psCurvenetP, psCurvenetE);
        psCurvenet->setColor({0.0f, 0.0f, 0.0f});   // Curvenet is grey
        psCurvenet->setMaterial("flat");
        psCurvenet->setTransparency(0.65);
        psCurvenet->setRadius(0.003);
        psCurvenet->setEnabled(true);

        editableCurvenet->convertControlsToPC(P);
        IO::convertVertsToMatrix(P, psControlsP);
        psControlsPC = polyscope::registerPointCloud("Controls", psControlsP);
        psControlsPC->setPointColor({0.9f, 0.2f, 0.1f});    // Controls are red
        psControlsPC->setMaterial("flat");
        psControlsPC->setPointRadius(0.02);
        psControlsCN->setEnabled(true);

        editableCurvenet->convertTangentsToPC(P);
        IO::convertVertsToMatrix(P, psTangentsP);
        psTangentsPC = polyscope::registerPointCloud("Tangents", psTangentsP);
        psTangentsPC->setPointColor({0.1f, 0.9f, 0.2f});    // Tangents are green
        psTangentsPC->setMaterial("flat");
        psTangentsPC->setPointRadius(0.012);
        psControlsCN->setEnabled(true);

        editableCurvenet->convertControlsAndTangentsToCN(P, psControlsAndTangentsE);
        IO::convertVertsToMatrix(P, psControlsAndTangents);
        psControlsCN = polyscope::registerCurveNetwork("Handles", psControlsAndTangents, psControlsAndTangentsE);
        psControlsCN->setColor({0.5f, 0.55f, 0.15f});   // Handles are yellow-ish
        psControlsCN->setMaterial("flat");
        psControlsCN->setTransparency(0.8);
        psControlsCN->setRadius(0.006);
        psControlsCN->setEnabled(true);
    
        // std::cout << "Curvenet updated with selection." << std::endl;
        recomputeCurvenet = false;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Too few arguments. Usage: ./profile_mover <input OBJ file path> <output file path>" << std::endl;
        return 1;
    }
    InputPath = argv[1];
    OutputPath = argv[2];

    // Initialize polyscope
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None; // Disable ground plane

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    polyscope::init();

    // Set camera view
    // polyscope::view::setUpDir(polyscope::UpDir::ZUp);      // Z up
    // polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront); // -Y forward

    // Set projection to orthographic
    polyscope::view::setProjectionMode(polyscope::ProjectionMode::Orthographic);

    // Load our mesh object
    std::cout << "\nLoading surface mesh file" << std::endl;
    IO::readOBJ(InputPath, V, F);
    IO::readOBJ(InputPath, psV, psF);

    // Register mesh with PS
    std::cout << "Registering Surface Mesh to Polyscope" << std::endl;
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", psV, psF);
    psMesh->setSurfaceColor({0.3f, 0.2f, 1.0f});

    // Create our own mesh object
    neutralMesh = std::make_unique<Mesh::mesh>(V, F);
    editableCurvenet = std::make_unique<Curvenet::curvenet>();
    // Alter PS edge permutation so we can match it
    // TODO: Check that this works.. PS doesn't seem to support non-triangular face edge picking
    // If not, fall back on computing the polyHE_t edge which is closest to the picked point
    // No need to clip t value in that case! Still, is quite costly...
    //std::vector<size_t> psEdgePerm = buildPolyscopeEdgePermutation(F, neutralMesh->getHEMesh());
    //psMesh->setEdgePermutation(psEdgePerm, neutralMesh->getHEMesh().numEdges());

    // Empty curvenet
    psCurvenet = polyscope::registerCurveNetwork("Curvenet", psCurvenetP, psCurvenetE);
    //psCurvenet->setEnabled(true);

    // Empty controls/tangent connectivity
    psControlsCN = polyscope::registerCurveNetwork("Handles", psControlsAndTangents, psControlsAndTangentsE);
    //psControlsCN->setEnabled(true);

    // Empty controls
    psControlsPC = polyscope::registerPointCloud("Controls", psControlsP);
    //psControlsPC->setEnabled(true);

    // Empty tangents
    psTangentsPC = polyscope::registerPointCloud("Tangents", psTangentsP);
    //psTangentsPC->setEnabled(true);

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}