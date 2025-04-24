
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/surface_centers.h"
#include "geometrycentral/surface/halfedge_mesh.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_mesh.h"
#include <iostream>
#include <cmath>
#include <random>
using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<HalfedgeMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
std::unique_ptr<VectorHeatMethodSolver> vhtSolver;
// so we can more easily pass these to different classes

ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

polyscope::SurfaceMesh* psMesh;
polyscope::VolumeMesh* cubeMesh;
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex

Vector<double> DELTA;

// Manage a list of sites for averages
struct SiteVert {
    Vertex vertex;
    float weight = 1.;
};
std::vector<SiteVert> siteVerts;
double vertexRadius,cubeLength;
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
void showGeodesicPath(std::vector<SurfacePoint> &path,int vidx){

    Vector3 xyzCoordinates{0,0,0};
    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    double tEdgevalue;
    for (auto point:path){
//
        xyzCoordinates=point.interpolate(geometry->inputVertexPositions);
        positions.push_back(xyzCoordinates);
        if (positions.size() >= 2) {
            edgeInds.push_back({positions.size() - 2, positions.size() - 1});
        }
    }
    polyscope::SurfaceGraphQuantity* geopathlines = psMesh->addSurfaceGraphQuantity(std::to_string(vidx), positions, edgeInds);
//    geopathlines->setEnabled(true);
    geopathlines->setRadius(vertexRadius*0.5);
    //  Random color for path line
    glm::vec3 RandomColor={randomReal(0,1), randomReal(0,1), randomReal(0,1)};
    geopathlines->setColor({0.57, 0.52, 0.75});

}
void addBoxMarker(double length){
    std::vector<Vector3> vertexPositions;
    std::vector<std::array<size_t, 8>> hexIndices;
    Vector3 basicvX,basicvY,basicvZ;
    int shift=0;
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        if (*it==polyscope::state::currVertexIndex) continue;
        Vector3 original=geometry->inputVertexPositions[*it];
        basicvX=geometry->vertexTangentBasis[*it][0];
        basicvY=geometry->vertexTangentBasis[*it][1];
        basicvZ=cross(basicvX,basicvY).normalize();
        //    add Cube vertices
        vertexPositions.push_back(original+length*(-basicvX-basicvY));    vertexPositions.push_back(original+length*(basicvX-basicvY));
        vertexPositions.push_back(original+length*(basicvX+basicvY));    vertexPositions.push_back(original+length*(-basicvX+basicvY));
        vertexPositions.push_back(original+length*(-basicvX-basicvY+basicvZ));    vertexPositions.push_back(original+length*(basicvX-basicvY+basicvZ));
        vertexPositions.push_back(original+length*(basicvX+basicvY+basicvZ));    vertexPositions.push_back(original+length*(-basicvX+basicvY+basicvZ));

        std::array<size_t, 8> hexcube;
        for (int idx=0;idx<8;++idx){
            hexcube[idx]=shift;
            shift+=1;
        }
        hexIndices.push_back(hexcube);

    }

    cubeMesh=polyscope::registerHexMesh("cubeMesh", vertexPositions, hexIndices);
    cubeMesh->setEnabled(true);
    cubeMesh->setColor({1.0, 0.45, 0.0});
}
void showSelected() {
    // Show selected vertices in yellow
    addBoxMarker(cubeLength);
    // Show the currently selected vertex in red.
    int currIdx = polyscope::state::currVertexIndex;
//    std::cout<<currIdx<<std::endl;
    if (currIdx != -1) {
        std::vector<Vector3> pos = {geometry->inputVertexPositions[currIdx]};
        currVert = psMesh->addSurfaceGraphQuantity("current vertex", pos, std::vector<std::array<size_t, 2>>());
        currVert->setEnabled(true);
        currVert->setRadius(vertexRadius);
        currVert->setColor({0.57, 0.52, 0.75});
    } else {
        currVert->setEnabled(false);
    }
}
void computeCenter() {

    if (vhtSolver == nullptr) {
        vhtSolver.reset(new VectorHeatMethodSolver(*geometry));
    }
    VertexData<double> dist(*mesh, 0.);
    for (std::set<size_t>::iterator it = polyscope::state::subset.vertices.begin();
         it != polyscope::state::subset.vertices.end(); ++it) {
        if (*it==polyscope::state::currVertexIndex) continue;
        dist[*it]=1;
    }
    int p=2;
//    SurfacePoint center = findCenter(*mesh,*geometry, *vhtSolver, dist, p);
//
//    // Visualize
//    Vector3 centerPos = center.interpolate(geometry->inputVertexPositions);
//    std::vector<Vector3> centerPosCloud{centerPos};
//    auto pointQ = polyscope::registerPointCloud("center", centerPosCloud);
//    pointQ->setPointRadius(vertexRadius);

//    Show path and points

//    std::vector<SurfacePoint> centerPosCloud;
    std::vector<Vector3> centerPosCloud;


    SurfacePoint initPoint=mesh->vertex(polyscope::state::currVertexIndex);
    std::vector<TraceGeodesicResult> searchPath=findCentervis(*mesh,*geometry, *vhtSolver, dist, p,initPoint);
    for (int pathIdx=0;pathIdx<(int)searchPath.size();++pathIdx){
        showGeodesicPath(searchPath[pathIdx].pathPoints,pathIdx);
//        centerPosCloud.push_back(searchPath[pathIdx].endPoint);
        centerPosCloud.push_back(searchPath[pathIdx].endPoint.interpolate(geometry->inputVertexPositions));
    }
    auto pointQ = polyscope::registerPointCloud("center", centerPosCloud);
    pointQ->setPointRadius(vertexRadius);

}
void redraw() {
    showSelected();
    psMesh->updateVertexPositions(geometry->inputVertexPositions);
    polyscope::requestRedraw();
}

void showExpMap(){
    if (vhtSolver == nullptr) {
        vhtSolver.reset(new VectorHeatMethodSolver(*geometry));
    }
    Vertex sourceV = mesh->vertex(polyscope::state::currVertexIndex);
    VertexData<Vector2> logmap= vhtSolver->computeLogMap(sourceV);
//    std::vector<std::array<double, 2>> V_uv(mesh->nVertices());
//    for (Vertex v : mesh->vertices()) {
//        V_uv[v.getIndex()] = {uvmap[v][0], uvmap[v][1]};
//    }

    auto checkerboard=psMesh->addLocalParameterizationQuantity("logmap", logmap);
    checkerboard->setEnabled(true);
    checkerboard->setStyle(polyscope::ParamVizStyle::LOCAL_CHECK);
    checkerboard->setCheckerColors(std::make_pair(glm::vec3{1.0, 0.45, 0.0}, glm::vec3{0.55, 0.27, 0.07}));
    checkerboard->setCheckerSize(0.05);

}
void functionCallback() {

    if (ImGui::Button("ExpontentialMap")) {
        if (polyscope::state::currVertexIndex!=-1) {
            showExpMap();
//            redraw();
        }
    }
    if (ImGui::Button("find center")) {
        computeCenter();
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
    psMesh = polyscope::registerSurfaceMesh("bunny", geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));



    // Initalize
    flipZ();
    double lengthScale = geometry->meanEdgeLength();
    vertexRadius = lengthScale * 0.2;
    cubeLength =lengthScale * 0.6;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({1.0, 1.0, 1.0}); // white
    currVert =
            psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    DELTA = Vector<double>::Zero(mesh->nVertices());
    geometry->requireVertexTangentBasis();
//    geometry->requireVertexNormals();


    polyscope::show();
//    delete mesh;
//    delete geometry;
//
    return EXIT_SUCCESS;

}