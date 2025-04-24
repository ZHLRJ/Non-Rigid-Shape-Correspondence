
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include <cmath>
#include <random>
using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
std::unique_ptr<VectorHeatMethodSolver> vhtSolver;
// so we can more easily pass these to different classes

ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
polyscope::SurfaceVertexParameterizationQuantity *checkerboard;
//polyscope::SurfaceVertexParameterizationQuantity
VertexData<Vector2> logmap;

double vertexRadius;
void flipZ() {
    // Rotate mesh 180 deg about up-axis on startup
    glm::mat4x4 rot = glm::rotate(glm::mat4x4(1.0f), static_cast<float>(PI), glm::vec3(0, 1, 0));
    for (Vertex v : mesh->vertices()) {

        Vector3 vec = geometry->inputVertexPositions[v];
        glm::vec4 rvec = {vec[0], vec[1], vec[2], 1.0};
        rvec = rot * rvec;
        geometry->inputVertexPositions[v] = {rvec[0], rvec[1], rvec[2]};
    }
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
}
void showSelected() {

    // Show the currently selected vertex in red.
    int currIdx = polyscope::state::currVertexIndex;
//    std::cout<<currIdx<<std::endl;
    if (currIdx != -1) {
        std::vector<Vector3> pos = {geometry->inputVertexPositions[currIdx]};
        currVert = psMesh->addSurfaceGraphQuantity("current vertex", pos, std::vector<std::array<size_t, 2>>());
        currVert->setEnabled(true);
        currVert->setRadius(vertexRadius);
        currVert->setColor({1.0, 0.0, 0.0});
    } else {
        currVert->setEnabled(false);
    }
}
void redraw() {
    showSelected();
    polyscope::requestRedraw();
}

void showExpMap(){
//    if (vhtSolver == nullptr) {
//        vhtSolver.reset(new VectorHeatMethodSolver(*geometry));
//    }
    logmap= vhtSolver->computeLogMap(mesh->vertex(polyscope::state::currVertexIndex));
//    std::vector<glm::vec2> V_uv(mesh->nVertices());
//    for (Vertex v : mesh->vertices()) {
//        V_uv[v.getIndex()] = glm::vec2{logmap[v][0], logmap[v][1]};
//    }
////    checkerboard=psMesh->addLocalParameterizationQuantity("logMap",logmap);
    checkerboard=psMesh->addLocalParameterizationQuantity("logMap",logmap);
    checkerboard->setEnabled(true);
    checkerboard->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    checkerboard->setCheckerColors(std::make_pair(glm::vec3{1.0, 0.45, 0.0}, glm::vec3{0.55, 0.27, 0.07}));
    checkerboard->setCheckerSize(0.05);

}


void functionCallback() {

    if (ImGui::Button("ExpontentialMap")) {
        if (polyscope::state::currVertexIndex!=-1) {
            showExpMap();

            }
        }

    if (ImGui::Button("Reset")) {

        polyscope::state::currVertexIndex = -1;
        psMesh->setSurfaceColor({1.0, 1.0, 1.0});
        redraw();
    }
}

int main(){
    //    std::complex<double> a(2,2);
    std::string filepath = "../../../../input/bunny.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);

    geometry=geometry_uptr.release();
    mesh=mesh_uptr.release();
    // Get indices for element picking
    polyscope::state::facePickIndStart = mesh->nVertices();
    polyscope::state::edgePickIndStart = polyscope::state::facePickIndStart + mesh->nFaces();
    polyscope::state::halfedgePickIndStart = polyscope::state::edgePickIndStart + mesh->nEdges();

    // construct a solver
    polyscope::init();
    // Set the callback function
    polyscope::state::userCallback = functionCallback;
    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh("logmap", geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));


    // Initalize
    flipZ();
    double lengthScale = geometry->meanEdgeLength();
    vhtSolver.reset(new VectorHeatMethodSolver(*geometry));

    vertexRadius = lengthScale * 0.2;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 1.0, 1.0}); // white
    currVert =
            psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());


    polyscope::show();
//    delete mesh;
//    delete geometry;
//
    return EXIT_SUCCESS;

}