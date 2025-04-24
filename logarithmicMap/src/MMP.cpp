// SIMPLICIAL COMPLEX

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

std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;

// == Geometry-central data
std::unique_ptr<HalfedgeMesh>* mesh2;
std::unique_ptr<ManifoldSurfaceMesh> *geometry2;

int main(){
    std::string filepath = "../../../../input/bunny.obj";
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
    geometry=geometry_uptr.release();
    mesh=mesh_uptr.release();

    VertexData<double> dist(*mesh, 0.);
    VectorHeatMethodSolver vhtSolver(*geometry);
//    SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
//                            const std::vector<Vertex>& vertexPts, int p)
//    SurfacePoint findCenter(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geom, VectorHeatMethodSolver& solver,
//                            const VertexData<double>& distribution, int p)

    SurfacePoint center = findCenter(*mesh,*geometry, vhtSolver, dist, 2);



//    Face f=mesh->face(1);
//
//    geometry->requireFaceTangentBasis();
//    Vector3 basicX=geometry->faceTangentBasis[f][0];
//    Vector3 basicY=geometry->faceTangentBasis[f][1];
//    std::cout<<basicX<<"  "<<basicY<<std::endl;
//    std::cout<<dot(basicY,basicX)<<std::endl;
//
//    geometry->requireVertexTangentBasis();
//
//    Vertex v=mesh->vertex(1);
//    Vector3 basicvX=geometry->vertexTangentBasis[v][0];
//    Vector3 basicvY=geometry->vertexTangentBasis[v][1];
//    std::cout<<basicvX<<"  "<<basicvY<<std::endl;
//    std::cout<<dot(basicvY,basicvX)<<std::endl;
//
//
//
}