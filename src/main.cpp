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
#include "profileformer/profileformer.hpp"
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

// Handles for discrete curvenet
std::vector<glm::vec3> psCurvenetP; // Point list
std::vector<std::array<size_t, 2>> psCurvenetE; // Edge List
polyscope::CurveNetwork* psCurvenet;

// Handles for controls (?)
// Not sure what we need visualization-wise
std::vector<glm::vec3> psControlsP; // Point list
polyscope::PointCloud* psControls;


// VARIABLES FOR PARSING AND WRITING FILES
std::string InputPath;

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



// ----------------- FUNCTIONS BEGIN HERE -------------------------

// Performs call to pre-computation of cut-mesh and operators
int computePrecomp() {
    return 1;
}

// Performs call to surface deformation and updates PS mesh
int computeDeformation() {
    return 1;
}


// Creates gizmo at vertex
void addGizmoAtVertex(int vert_idx) {
    glm::vec3 startpos = psV[vert_idx];

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
    psV[vert_idx] = new_pos;
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
    Utils::copyPositions(T, psT);
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", psV, psT);
    return;
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
        polyscope::VolumeMeshPickResult meshPick = psMesh->interpretPickResult(pick);

        // TODO: We should be able to click the mesh arbitrarily, not just at vertices
        if (meshPick.elementType == polyscope::VolumeMeshElement::VERTEX) {
            selectedVertex = static_cast<int>(meshPick.index);
        } else {
            std::cout << "Did not click on mesh." << std::endl;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Too few arguments. Usage: ./profile_former <OBJ file path>" << std::endl;
        return 1;
    }
    InputPath = argv[1];

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
    igl::readOBJ(InputPath, V, T);

    // Register tet mesh with PS
    std::cout << "Registering Surface Mesh to Polyscope" << std::endl;
    psMesh = polyscope::registerSurfaceMesh("Surface Mesh", V, T);

    // Empty curvenet
    psCurvenet = polyscope::registerCurveNetwork("Discrete Curvenet", psCurvenetP, psCurvenetE);
    //psCurvenet->setEnabled(true);

    // Empty controls
    psControls = polyscope::registerPointCloud("Controls", psCurvenetP);
    //psControls->setEnabled(true);

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}