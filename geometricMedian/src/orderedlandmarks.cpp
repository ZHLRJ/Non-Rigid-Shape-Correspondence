//
// Created by mars_zhang on 10/18/23.
//

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/surface_centers.h"
#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"

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
std::unique_ptr<HeatMethodDistanceSolver> heatSolver;
// so we can more easily pass these to different classes

ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

polyscope::SurfaceMesh* psMesh;
polyscope::VolumeMesh* cubeMesh;
polyscope::PointCloud* PCloud;
polyscope::PointCloud* salientCloud;
polyscope::PointCloud* pointQ;
polyscope::SurfaceGraphQuantity* currVert; // currently active vertex

//Some global variable
VertexData<double> dist;
std::vector<Vector3> CenterPos;
std::vector<SurfacePoint> CenterSurfPoints;
std::vector<std::pair<double,size_t>> geoDistance;
SurfacePoint centerOfKarcher;
double maxDualArea=0;

// Manage a list of sites for averages
struct SiteVert {
    Vertex vertex;
    float weight = 1.;
};
std::vector<SiteVert> siteVerts;
double vertexRadius,cubeLength;

int nRandVerts=2000;
int pickedSalient=5;
int numKarcherMean =10;
std::vector<std::pair<double,std::vector<double>>> salientInfo;
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
    geopathlines->setRadius(vertexRadius*0.3);
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
std::vector<int> shuffleVector(size_t N,bool fixedSeed=false){
    //    randome pick N points
    int seed = 0;
    // set some values:
    std::vector<int> shuffleVector(N);
    std::iota(shuffleVector.begin(),shuffleVector.end(),0);

    if (fixedSeed){
        std::mt19937 g(seed);
        std::shuffle(shuffleVector.begin(), shuffleVector.end(), g);
    }
    else {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(shuffleVector.begin(), shuffleVector.end(), g);
    }
    return shuffleVector;
}
void randomSelected(){
    std::vector<Vector3> vertPos;
    std::vector<int> v = shuffleVector(mesh->nVertices());
    nRandVerts = std::min(nRandVerts,(int)mesh->nVertices());

//    std::vector<std::pair<size_t, int>> vCount;
//    for (int iV = 0; iV < nRandVerts; iV++) {
//        vCount.push_back(std::make_pair(v[iV], 2));
//        }
//    psMesh->addVertexCountQuantity("Random selected Vertices",vCount);
    // Show the picked vertex.
    dist=VertexData<double> (*mesh, 0.);
    for (auto it=0; it <nRandVerts; ++it) {
        vertPos.push_back(geometry->inputVertexPositions[v[it]]);
//        dist[mesh->vertex(v[it])] = 1;

//        dist[mesh->vertex(v[it])] = geometry->vertexDualArea(mesh->vertex(v[it]))/maxDualArea;
        dist[mesh->vertex(v[it])] = geometry->vertexDualArea(mesh->vertex(v[it]));

//        printf("Vertex %d area: %.4f \n",v[it],dist[mesh->vertex(v[it])]);

    }
//    //    use pointCloud structure to visualize picked vertices
    PCloud = polyscope::registerPointCloud("random points", vertPos);
    PCloud->setEnabled(true);
    PCloud->setPointRadius(vertexRadius*0.1);
    PCloud->setPointColor({1.0, 0.65, 0.0});

    psMesh->addVertexScalarQuantity("density map",dist);

}
void showSelected() {
    // Show selected vertices in yellow
//    addBoxMarker(cubeLength);
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
//    VertexData<double> dist(*mesh,1.);
//    for (Vertex v:mesh->vertices()) {
//        int vid=v.getIndex();
//        dist[vid]=1+geometry->inputVertexPositions[vid].x;
//    }
    int p=2;
    CenterPos.clear();
    CenterSurfPoints.clear();
    VertexData<double> meanDist=VertexData<double> (*mesh, 0.);

    for (int i=0;i<numKarcherMean;i++)
    {
        SurfacePoint center = findCenter(*mesh, *geometry, *vhtSolver, dist, p);
        CenterPos.push_back(center.interpolate(geometry->inputVertexPositions)) ;
        CenterSurfPoints.push_back(center) ;
//        std::cout<<"The nearest vertex "<<center.nearestVertex()<<std::endl;
        meanDist[center.nearestVertex()] = 1;

    }
//    currVert = psMesh->addSurfaceGraphQuantity("TESTKarcherMean", CenterPos, std::vector<std::array<size_t, 2>>());
//    currVert->setEnabled(true);
//    currVert->setRadius(vertexRadius);
//    currVert->setColor({0.57, 0.52, 0.75});
//    currVert->setColor({0.435, 0.482, 0.631});

//    Mean of all Karcher means
    centerOfKarcher = findCenter(*mesh, *geometry, *vhtSolver, meanDist, p);


    polyscope::SurfaceGraphQuantity * RefineCenter =
            psMesh->addSurfaceGraphQuantity("RefineCenter",std::vector<Vector3> {centerOfKarcher.interpolate(geometry->inputVertexPositions)}, std::vector<std::array<size_t, 2>>());

    RefineCenter->setEnabled(true);
    RefineCenter->setRadius(vertexRadius);
//    currVert->setColor({0.57, 0.52, 0.75});
    RefineCenter->setColor({1, 0.482, 0.631});
}

void refineGeoMedian(){
    if (vhtSolver == nullptr) {
        vhtSolver.reset(new VectorHeatMethodSolver(*geometry));
    }

    VertexData<Vector2> logmap = vhtSolver->computeLogMap(centerOfKarcher);
    std::vector<std::pair<double,int>> distance;
//    distance.resize(CenterPos.size());

    for (int it=0; it <CenterPos.size(); ++it) {
        Vertex v= CenterSurfPoints[it].nearestVertex();

        Vector2 pointCoord = logmap[v];
        double dist2 = pointCoord.norm2();
        distance.push_back(std::pair<double,int> (dist2,it));
    }
    std::sort(distance.begin(),distance.end());

    std::vector<Vector3> RefineCenterPos;
//    CenterPos.size()/2
    for (int it=0; it <4; ++it){
        RefineCenterPos.push_back(CenterPos[distance[it].second]);
//        printf(" ( %f , %d ) \n",distance[it].first,distance[it].second);
    }

    CenterPos = RefineCenterPos;
    currVert = psMesh->addSurfaceGraphQuantity("KarcherMean", CenterPos, std::vector<std::array<size_t, 2>>());
    currVert->setEnabled(true);
    currVert->setRadius(vertexRadius);
//    currVert->setColor({0.57, 0.52, 0.75});
    currVert->setColor({0.435, 0.482, 0.631});

//    for (int it=0; it <CenterPos.size(); ++it){
//        printf(" ( %f , %d ) \n",distance[it].first,distance[it].second);
//    }


}
void findGreatestDistance(){
    if (vhtSolver == nullptr) {
        vhtSolver.reset(new VectorHeatMethodSolver(*geometry));
    }
    geoDistance.resize(mesh->nVertices());
    std::vector<std::pair<double,size_t>>::iterator item ;
    for (item= geoDistance.begin();item<geoDistance.end();++item){
        item->first = 0;
        item->second = item -geoDistance.begin();
    }
//    std::fill(geoDistance.begin(), geoDistance.end(), 0);
//    distToSource.resize(mesh->nVertices());
    std::cout <<CenterPos.size()<<std::endl;
    for (size_t idx=0;idx<CenterPos.size();++idx) {
//        std::cout<<idx<<CenterPos[idx]<<std::endl;
        VertexData<Vector2> logmap = vhtSolver->computeLogMap(CenterSurfPoints[idx]);
        for (Vertex v : mesh->vertices()) {

            Vector2 pointCoord = logmap[v];
            double dist2 = pointCoord.norm2();
            geoDistance[v.getIndex()].first+=dist2;
//            distToSource[v.getIndex()] = dist2;

        }
//        std::vector<double>::iterator largest = std::max_element(distToSource.begin(),distToSource.end());
//        printf("largest dis %.2f \n",*largest);

    }
    std::sort(geoDistance.rbegin(),geoDistance.rend());
//    test
//    std::vector<std::pair<double,size_t>> test;
//    test.assign(geoDistance.begin(),geoDistance.begin()+5);
//    for (auto item:test){
//        printf("The index %zd; distance %.2f\n",item.second,item.first);
//
//    }

}
std::vector<double> calcdistToSource(int vidx){
    std::vector<double> distToSource(mesh->nVertices(),0);
    VertexData<Vector2> logmap = vhtSolver->computeLogMap(mesh->vertex(vidx));

    for (Vertex v : mesh->vertices()) {
        Vector2 pointCoord = logmap[v];
        double dist2 = pointCoord.norm2();
        distToSource[v.getIndex()] = dist2;
    }
//    std::vector<double>::iterator largest = std::max_element(distToSource.begin(),distToSource.end());
//    printf("largest dis %.2f \n",*largest);
    return distToSource;
}
bool distancetest(int vidx){
 for (auto item:salientInfo){
     if (0.8*item.first>item.second[vidx])
         return false;
 }
 return true;
}
void showSalient(int picked =1){
    std::vector<Vector3> salientPoints;
    std::vector<Vector3> unstaisfysalientPoints;
//    if (heatSolver == nullptr) {
//        heatSolver.reset(new HeatMethodDistanceSolver(*geometry));
//    }
    int count =0;
//    for (int i=0; i<picked;i++)
//        salientPoints.push_back(geometry->inputVertexPositions[geoDistance[i].second]);

    salientInfo.clear();
    for (int i=0; i<mesh->nVertices();++i) {
        double dis = geoDistance[i].first;
        size_t vidx = geoDistance[i].second;
        if (salientPoints.empty())
        {
//            std::cout<<"Empty"<<std::endl;

            std::vector<double> distToSource = calcdistToSource(vidx);
            salientPoints.push_back(geometry->inputVertexPositions[vidx]);
            salientInfo.push_back(std::make_pair(dis/CenterPos.size(),distToSource));
            count++;
            printf("Add the first vertex as salient %zd\n",vidx);

        }
        else if (distancetest(vidx) ) {

            std::vector<double> distToSource = calcdistToSource(vidx);
            salientPoints.push_back(geometry->inputVertexPositions[vidx]);
            salientInfo.push_back(std::make_pair(dis/CenterPos.size(), distToSource));
            count++;
            printf("Add the vertex as salient %zd\n",vidx);

        }
        else {
            unstaisfysalientPoints.push_back(geometry->inputVertexPositions[vidx]);
        }

        if (count>=picked)
            break;
    }


    std::vector<double> colorScalar(salientPoints.size());
    for (size_t i = 0; i < salientPoints.size(); i++) {
        colorScalar[i] = i+2;
    }
    salientCloud = polyscope::registerPointCloud("salient Points", salientPoints);
    salientCloud->setEnabled(true);
    salientCloud->setPointRadius(vertexRadius*2);
//    salientCloud->addColorQuantity("order color",);
    salientCloud->addScalarQuantity("salient color order", colorScalar);

    currVert = psMesh->addSurfaceGraphQuantity("unstaisfysalientPoints", unstaisfysalientPoints, std::vector<std::array<size_t, 2>>());
    currVert->setEnabled(false);
    currVert->setRadius(vertexRadius);
    currVert->setColor({0.57, 0.0, 0.75});

}
void functionCallback() {
    //some variable
    ImGui::InputInt("N verts",&nRandVerts);
    ImGui::InputInt("N salient",&pickedSalient);

    if (ImGui::Button("random select")){
        randomSelected();
    }
    if (ImGui::Button("find GeoMedian")) {
        auto start = std::chrono::steady_clock::now();
        calculateGeoMedian();
        auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_seconds{end - start};
        std::cout << "calculateGeoMedian() use " << elapsed_seconds.count() << " s" << std::endl;

    }

    if (ImGui::Button("refineGeoMedian")) {
        auto start = std::chrono::steady_clock::now();
        refineGeoMedian();
        auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_seconds{end - start};
        std::cout << "refineGeoMedian() use " << elapsed_seconds.count() << " s" << std::endl;

    }
    if (ImGui::Button("find Greatest Distance")) {
        auto start = std::chrono::steady_clock::now();
        findGreatestDistance();
        auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_seconds{end - start};
        std::cout << "findGreatestDistance() use " << elapsed_seconds.count() << " s" << std::endl;

    }
    if (ImGui::Button("show salient Points")) {
        auto start = std::chrono::steady_clock::now();

        showSalient(pickedSalient);
        auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_seconds{end - start};
        std::cout << "showSalient(pickedSalient) use " << elapsed_seconds.count() << " s" << std::endl;

    }



//    if (ImGui::Button("ExpontentialMap")) {
//        if (polyscope::state::currVertexIndex!=-1) {
//            showExpMap();
////            redraw();
//        }
//    }
//    if (ImGui::Button("find center")) {
//        computeCenter();
//    }
//    if (ImGui::Button("find GeoMedian")) {
//        calculateGeoMedian();
//    }
//
//    if (ImGui::Button("Reset")) {
//        polyscope::state::currVertexIndex = -1;
////        psMesh->setSurfaceColor({1.0, 1.0, 1.0});

//        redraw();
//    }
}

int main(){

//    std::string filepath="/Users/mars_zhang/Downloads/Github/shapeMatching/Models/MPI-FAUST/training/registrations/tr_reg_005.ply";
    std::string filepath="../../../../input/vtx_5k/cat3.off";
//    std::string filepath= "/Users/mars_zhang/Downloads/ZHL_Paper/IEEEVIS/code/Iandmarks/data/SGA18_orientation_BCICP_dataset-master/Dataset/FAUST/vtx_5k/tr_reg_002.off";
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
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::groundPlaneHeightFactor = 0.02; // adjust the plane height
    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh("FAUST", geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));



    // Initalize
    geometry->requireVertexDualAreas();
    VertexData<double> dualArea = geometry->vertexDualAreas;
    for( Vertex v : mesh->vertices() ){
        if (dualArea[v]>maxDualArea) {
            maxDualArea = dualArea[v];
        }
    }
    printf("The maximum vertex dual Area : %.4f\n",maxDualArea);

    flipZ();
    double lengthScale = geometry->meanEdgeLength();  // 0.016
    printf("The average lengthScale %.3f",lengthScale);
    vertexRadius = 0.016 ;
    cubeLength =lengthScale * 0.6;
    psMesh->setSmoothShade(true);
    psMesh->setSurfaceColor({0.95, 0.94, 0.94}); // white
    currVert =
            psMesh->addSurfaceGraphQuantity("current vertex", std::vector<Vector3>(), std::vector<std::array<size_t, 2>>());
    geometry->requireVertexTangentBasis();
//    geometry->requireVertexNormals(); polyscope::view::upDir = UpDir::ZUp;

    polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
    polyscope::show();
//    delete mesh;
//    delete geometry;
//
    return EXIT_SUCCESS;

}
