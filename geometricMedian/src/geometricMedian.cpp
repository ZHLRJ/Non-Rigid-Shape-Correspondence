//
// Created by mars_zhang on 2/22/23.
//

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
polyscope::PointCloud* PCloud;
polyscope::PointCloud* pointQ;
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex
VertexData<double> dist;

// Manage a list of sites for averages
struct SiteVert {
    Vertex vertex;
    float weight = 1.;
};
std::vector<SiteVert> siteVerts;
double vertexRadius,cubeLength;
int nRandVerts=100;
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
void randomSelected(bool randomseed=false){
    //    randome pick N points
    std::random_device rd;
    std::mt19937 g(rd());
    if (randomseed) {
        int seed = 0;
        std::mt19937 g(seed);
    }
    std::vector<int> v (mesh->nVertices());
    std::iota(v.begin(),v.end(),0);
    std::shuffle(v.begin(), v.end(), g);

    // Show the picked vertex.
    std::vector<Vector3> vertPos;
    std::vector<std::array<size_t, 2>> vertInd;
    for (size_t it=0; it <nRandVerts; ++it) {
        vertPos.push_back(geometry->inputVertexPositions[v[it]]);
        dist[mesh->vertex(it)]=1;
    }
//    polyscope::SurfaceGraphQuantity* showVerts = psMesh->addSurfaceGraphQuantity("random picked vertices", vertPos, vertInd);
//    showVerts->setEnabled(true);
//    showVerts->setRadius(vertexRadius);
//    showVerts->setColor({1.0, 0.65, 0.0});
    //    use pointCloud structure to save
    PCloud = polyscope::registerPointCloud("random points", vertPos);
    PCloud->setEnabled(true);
    PCloud->setPointRadius(vertexRadius);
    PCloud->setPointColor({1.0, 0.65, 0.0});

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
//        dist[*it]=geometry->inputVertexPositions[*it].x;
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
    std::vector<TraceGeodesicResult> searchPath=findCenterwithPath(*mesh,*geometry, *vhtSolver, dist, p,initPoint);
    for (int pathIdx=0;pathIdx<(int)searchPath.size();++pathIdx){
        showGeodesicPath(searchPath[pathIdx].pathPoints,pathIdx);
//        centerPosCloud.push_back(searchPath[pathIdx].endPoint);
        centerPosCloud.push_back(searchPath[pathIdx].endPoint.interpolate(geometry->inputVertexPositions));
    }
    pointQ = polyscope::registerPointCloud("center", centerPosCloud);
    pointQ->setPointRadius(vertexRadius);
}
void redraw() {
    PCloud->remove();
    pointQ->remove();
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
void calculateGeoMedian(){
    if (vhtSolver == nullptr) {
        vhtSolver.reset(new VectorHeatMethodSolver(*geometry));
    }
    //  Geometric Median
    VertexData<double> dist(*mesh,1.);
//    for (Vertex v:mesh->vertices()) {
//        int vid=v.getIndex();
//        dist[vid]=1+geometry->inputVertexPositions[vid].x;
//    }
    int p=2;
    if (p==2){
        std::vector<Vector3> pos;
        for (int i=0;i<15;i++) {
            SurfacePoint center = findCenter(*mesh, *geometry, *vhtSolver, dist, p);
             pos.push_back(center.interpolate(geometry->inputVertexPositions)) ;
        }
        currVert = psMesh->addSurfaceGraphQuantity("current vertex", pos, std::vector<std::array<size_t, 2>>());
        currVert->setEnabled(true);
        currVert->setRadius(vertexRadius);
        currVert->setColor({0.57, 0.52, 0.75});
    }



//    std::vector<Vector3> centerPosCloud;
//    if (dist.toVector().maxCoeff()!=0){
//        SurfacePoint initPoint=mesh->vertex(polyscope::state::currVertexIndex);
////        SurfacePoint initPoint=mesh->vertex(1050);
//
//        std::vector<TraceGeodesicResult> searchPath=findCenterwithPath(*mesh,*geometry, *vhtSolver, dist, p,initPoint);
//
//        for (int pathIdx=0;pathIdx<(int)searchPath.size();++pathIdx){
//            showGeodesicPath(searchPath[pathIdx].pathPoints,pathIdx);
//    //        centerPosCloud.push_back(searchPath[pathIdx].endPoint);
//            centerPosCloud.push_back(searchPath[pathIdx].endPoint.interpolate(geometry->inputVertexPositions));
//        }
//        pointQ = polyscope::registerPointCloud("Geocenter", centerPosCloud);
//        pointQ->setPointRadius(vertexRadius);
//    }
}
void functionCallback() {
    //some variable
    ImGui::InputInt("N verts",&nRandVerts);
//    if (ImGui::Button("ExpontentialMap")) {
//        if (polyscope::state::currVertexIndex!=-1) {
//            showExpMap();
////            redraw();
//        }
//    }
//    if (ImGui::Button("find center")) {
//        computeCenter();
//    }
    if (ImGui::Button("find GeoMedian")) {
        calculateGeoMedian();
    }

    if (ImGui::Button("random select")){
        randomSelected(true);
    }

    if (ImGui::Button("Reset")) {
        polyscope::state::currVertexIndex = -1;
//        psMesh->setSurfaceColor({1.0, 1.0, 1.0});

        redraw();
    }
}

int main(){
    //    std::complex<double> a(2,2);
//    std::string filepath = "../../../../input/tr_reg_002.obj";
    std::string filepath="../../../../input/FAUST_test/tr_reg_080.ply";

//     std::string filepath = "../../../../input/smallDataset/centaur3.obj";
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
    geometry->requireVertexTangentBasis();
//    geometry->requireVertexNormals();
    dist=VertexData<double> (*mesh, 0.);


    polyscope::show();
//    delete mesh;
//    delete geometry;
//
    return EXIT_SUCCESS;

}
