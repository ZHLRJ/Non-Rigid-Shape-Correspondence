#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_centers.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <sstream>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<HalfedgeMesh> mesh;
ManifoldSurfaceMesh *geometry;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

// Some algorithm parameters
float tCoef = 1.0;
std::unique_ptr<VectorHeatMethodSolver> solver;
int vertexInd = 0;
int pCenter = 2;

// Manage a list of sources
struct SourceVert {
  Vertex vertex;
  float scalarVal = 1.;
  float vectorMag = 1.;
  float vectorAngleRad = 0.;
};
std::vector<SourceVert> sourcePoints;

// Manage a list of sites for averages
struct SiteVert {
  Vertex vertex;
  float weight = 1.;
};
std::vector<SiteVert> siteVerts;

bool vizFirstRun = true;
void updateSourceSetViz() {

  // Scalar balls around sources
  std::vector<std::pair<size_t, double>> sourcePairs;
  for (SourceVert& s : sourcePoints) {
    size_t ind = geometry->vertexIndices[s.vertex];
    sourcePairs.emplace_back(ind, s.scalarVal);
  }
  auto scalarQ = polyscope::getSurfaceMesh()->addVertexIsolatedScalarQuantity("source scalars", sourcePairs);
  scalarQ->pointRadius *= 2.;
  scalarQ->cMap = "reds";
  if (vizFirstRun) {
    scalarQ->setEnabled(true);
  }

  // Vectors at sources
  VertexData<Vector2> sourceVectors(*mesh, Vector2::zero());
  for (SourceVert& s : sourcePoints) {
    sourceVectors[s.vertex] = Vector2::fromAngle(s.vectorAngleRad) * s.vectorMag;
  }
  auto vectorQ = polyscope::getSurfaceMesh()->addVertexIntrinsicVectorQuantity("source vectors", sourceVectors);
  vectorQ->setVectorLengthScale(.05);
  vectorQ->setVectorRadius(.005);
  vectorQ->setVectorColor(glm::vec3{227 / 255., 52 / 255., 28 / 255.});
  if (vizFirstRun) {
    vectorQ->setEnabled(true);
  }

  vizFirstRun = false;
}

bool vizFirstRunSite = true;
void updateSiteSetViz() {

  // Scalar balls around sources
  std::vector<std::pair<size_t, double>> sourcePairs;
  for (SiteVert& s : siteVerts) {
    size_t ind = geometry->vertexIndices[s.vertex];
    sourcePairs.emplace_back(ind, s.weight);
  }
  auto scalarQ = polyscope::getSurfaceMesh()->addVertexIsolatedScalarQuantity("averaging sites", sourcePairs);
  scalarQ->pointRadius *= 2.;
  scalarQ->cMap = "blues";
  if (vizFirstRunSite) {
    scalarQ->setEnabled(true);
  }

  vizFirstRunSite = false;
}
void addVertexSource(size_t ind) {
  Vertex v = mesh->vertex(ind);

  // Make sure not already used
  for (SourceVert& s : sourcePoints) {
    if (s.vertex == v) {
      std::stringstream ss;
      ss << "Vertex " << v;
      std::string vStr = ss.str();
      polyscope::warning("Vertex " + vStr + " is already a source");
      return;
    }
  }

  SourceVert newV;
  newV.vertex = v;
  sourcePoints.push_back(newV);
  updateSourceSetViz();
}
void addVertexSite(size_t ind) {
  Vertex v = mesh->vertex(ind);

  // Make sure not already used
  for (SiteVert& s : siteVerts) {
    if (s.vertex == v) {
      std::stringstream ss;
      ss << "Vertex " << v;
      std::string vStr = ss.str();
      polyscope::warning("Vertex " + vStr + " is already a site");
      return;
    }
  }

  SiteVert newV;
  newV.vertex = v;
  newV.weight = 1.0;
  siteVerts.push_back(newV);
  updateSiteSetViz();
}

void scalarExtension() {
  if (solver == nullptr) {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }

  if (sourcePoints.size() == 0) {
    polyscope::warning("no source points set");
    return;
  }

  std::vector<std::tuple<SurfacePoint, double>> points;
  for (SourceVert& s : sourcePoints) {
    points.emplace_back(s.vertex, s.scalarVal);
  }

  VertexData<double> scalarExtension = solver->extendScalar(points);

  auto psScalar = psMesh->addVertexScalarQuantity("scalar extension", scalarExtension);
  psScalar->setEnabled(true);
}

void vectorTransport() {
  if (solver == nullptr) {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }

  if (sourcePoints.size() == 0) {
    polyscope::warning("no source points set");
    return;
  }

  std::vector<std::tuple<SurfacePoint, Vector2>> points;
  for (SourceVert& s : sourcePoints) {
    points.emplace_back(s.vertex, Vector2::fromAngle(s.vectorAngleRad) * s.vectorMag);
  }
  VertexData<Vector2> vectorExtension = solver->transportTangentVectors(points);

  auto psVec = psMesh->addVertexIntrinsicVectorQuantity("vector extension", vectorExtension);
  psVec->setEnabled(true);
}

void computeLogMap() {
  if (solver == nullptr) {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }
  if (sourcePoints.size() == 0) {
    polyscope::warning("must select a source");
    return;
  }

  Vertex sourceV = sourcePoints[0].vertex;
  VertexData<Vector2> logmap = solver->computeLogMap(sourceV);

  auto psLogmap = psMesh->addLocalParameterizationQuantity("logmap", logmap);
  psLogmap->setEnabled(true);
}


void computeCenter() {
  if (!(pCenter == 1 || pCenter == 2)) {
    polyscope::warning("p must be 1 or 2");
    return;
  }
  if (siteVerts.size() == 0) {
    polyscope::warning("must select at least one site");
    return;
  }

  if (solver == nullptr) {
    solver.reset(new VectorHeatMethodSolver(*geometry, tCoef));
  }

  // Build distribution
  VertexData<double> dist(*mesh, 0.);
  for (SiteVert& s : siteVerts) {
    dist[s.vertex] += s.weight;
  }

  SurfacePoint center = findCenter(*geometry, *solver, dist, pCenter);

  // Visualize
  Vector3 centerPos = center.interpolate(geometry->inputVertexPositions);
  std::vector<Vector3> centerPosCloud{centerPos};
  auto pointQ = polyscope::registerPointCloud("center", centerPosCloud);
  pointQ->setPointRadius(5.);
}

void buildPointsMenu() {

  bool anyChanged = false;

  ImGui::PushItemWidth(200);

  int id = 0;
  int eraseInd = -1;
  for (SourceVert& s : sourcePoints) {
    std::stringstream ss;
    ss << "Vertex " << s.vertex;
    std::string vStr = ss.str();
    ImGui::PushID(vStr.c_str());

    ImGui::TextUnformatted(vStr.c_str());

    ImGui::SameLine();
    if (ImGui::Button("delete")) {
      eraseInd = id;
      anyChanged = true;
    }
    ImGui::Indent();

    if (ImGui::InputFloat("scalar value", &s.scalarVal)) anyChanged = true;
    if (ImGui::InputFloat("vector mag", &s.vectorMag)) anyChanged = true;
    if (ImGui::SliderAngle("vector angle", &s.vectorAngleRad)) anyChanged = true;

    ImGui::Unindent();
    ImGui::PopID();
  }
  ImGui::PopItemWidth();

  // actually do erase, if requested
  if (eraseInd != -1) {
    sourcePoints.erase(sourcePoints.begin() + eraseInd);
  }

  if (ImGui::Button("add point")) {
    long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
    if (pickVert >= 0) {
      addVertexSource(pickVert);
      anyChanged = true;
    }
  }

  if (anyChanged) {
    updateSourceSetViz();
  }
}

void buildSitesMenu() {

  bool anyChanged = false;

  ImGui::PushItemWidth(200);

  int id = 0;
  int eraseInd = -1;
  for (SiteVert& s : siteVerts) {
    std::stringstream ss;
    ss << "Vertex " << s.vertex;
    std::string vStr = ss.str();
    ImGui::PushID(vStr.c_str());

    ImGui::TextUnformatted(vStr.c_str());

    ImGui::SameLine();
    if (ImGui::Button("delete")) {
      eraseInd = id;
      anyChanged = true;
    }
    ImGui::Indent();

    if (ImGui::InputFloat("weight", &s.weight)) anyChanged = true;

    ImGui::Unindent();
    ImGui::PopID();
  }
  ImGui::PopItemWidth();

  // actually do erase, if requested
  if (eraseInd != -1) {
    siteVerts.erase(siteVerts.begin() + eraseInd);
  }

  if (ImGui::Button("add site")) {
    long long int pickVert = polyscope::getSurfaceMesh()->selectVertex();
    if (pickVert >= 0) {
      addVertexSite(pickVert);
      anyChanged = true;
    }
  }

  if (anyChanged) {
    updateSiteSetViz();
  }
}

void myCallback() {

  ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;
  if (ImGui::BeginTabBar("MyTabBar", tab_bar_flags)) {
    if (ImGui::BeginTabItem("Basic algorithm")) {

      ImGui::TextUnformatted("Algorithm options:");
      ImGui::PushItemWidth(100);
      if (ImGui::InputFloat("tCoef", &tCoef)) {
        solver.reset();
      }
      ImGui::PopItemWidth();

      // Build the list of source points
      if (ImGui::TreeNode("select source points")) {
        buildPointsMenu();
        ImGui::TreePop();
      }

      if (ImGui::Button("run scalar extension")) {
        scalarExtension();
      }

      if (ImGui::Button("run vector transport")) {
        vectorTransport();
      }

      if (ImGui::Button("compute log map")) {
        computeLogMap();
      }


      ImGui::EndTabItem();
    }
    if (ImGui::BeginTabItem("Centers")) {

      if (ImGui::TreeNode("select sites to compute center of")) {
        buildSitesMenu();
        ImGui::TreePop();
      }


      ImGui::PushItemWidth(200);
      ImGui::InputInt("p norm", &pCenter);
      ImGui::PopItemWidth();

      if (ImGui::Button("find center")) {
        computeCenter();
      }

      ImGui::EndTabItem();
    }
    ImGui::EndTabBar();
  }
}

int main(int argc, char** argv) {

  // Configure the argument parser
  args::ArgumentParser parser("A demo of the Vector Heat Method");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help&) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load mesh
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(args::get(inputFilename)),
                                          geometry->inputVertexPositions, mesh->getFaceVertexList(),
                                          polyscopePermutations(*mesh));


  // Set vertex tangent spaces
  geometry->requireVertexTangentBasis();
  VertexData<Vector3> vBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  polyscope::getSurfaceMesh()->setVertexTangentBasisX(vBasisX);

  // Set face tangent spaces
  geometry->requireFaceTangentBasis();
  FaceData<Vector3> fBasisX(*mesh);
  for (Face f : mesh->faces()) {
    fBasisX[f] = geometry->faceTangentBasis[f][0];
  }
  polyscope::getSurfaceMesh()->setFaceTangentBasisX(fBasisX);

  // To start, pick two vertices as sources
  geometry->requireVertexIndices();
  addVertexSource(0);
  addVertexSource(mesh->nVertices() / 2);
  sourcePoints[1].scalarVal = 3.0;


  addVertexSite(0);
  addVertexSite(2 * mesh->nVertices() / 3);
  addVertexSite(mesh->nVertices() / 3);


  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
