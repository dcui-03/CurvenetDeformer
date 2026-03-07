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

    const auto controlsD = Utils::buildControlCloud(PF->getNeutralCurvenet());
    psControlsP.clear();
    psControlsP.reserve(controlsD.size());
    for (const auto& p : controlsD) {
        psControlsP.emplace_back(static_cast<float>(p[0]),
                                 static_cast<float>(p[1]),
                                 static_cast<float>(p[2]));
    }

    psCurvenet = polyscope::registerCurveNetwork("Discrete Curvenet", psCurvenetP, psCurvenetE);
    psControls = polyscope::registerPointCloud("Controls", psControlsP);
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
        polyscope::SurfaceMeshPickResult meshPick = psMesh->interpretPickResult(pick);

        // TODO: We should be able to click the mesh arbitrarily, not just at vertices
        if (meshPick.elementType == polyscope::MeshElement::VERTEX) {
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
    PF = std::make_unique<ProfileFormer::profileformer>(V, T, cnInput.controlP, cnInput.curveC);
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
}
