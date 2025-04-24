//
//MMP method
//

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include <iostream>
#include <cmath>
#include <random>
using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh>  mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::vector<Vector3> KN_vectors;

int main(){
//    std::complex<double> a(2,2);
//    std::string filepath = "../../../../input/bunny.obj";
//    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);
//    VertexData<Vector2> VectorHeatSolver::computeLogMap(const Vertex& sourceVert)

}