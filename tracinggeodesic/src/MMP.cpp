//
// Created by mars_zhang on 1/12/23.
//
#include "geometrycentral/surface/exact_geodesics.h"
#include "geometrycentral/surface/meshio.h"
#include <iostream>
using namespace geometrycentral;
using namespace geometrycentral::surface;
// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh_uptr;
std::unique_ptr<VertexPositionGeometry> geometry_uptr;
//// so we can more easily pass these to different classes
ManifoldSurfaceMesh* mesh;
VertexPositionGeometry* geometry;
int main(){
    std::string filename = "../../../../input/bunny.obj";
    // Load mesh
    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filename);
    mesh = mesh_uptr.release();
    geometry = geometry_uptr.release();

    // Create the GeodesicAlgorithmExact
    GeodesicAlgorithmExact mmp(*mesh, *geometry);

    // Pick a few points as the source set
    std::vector<SurfacePoint> sourcePoints;
    Vertex v = mesh->vertex(0);
    sourcePoints.push_back(SurfacePoint(v));

    // Run MMP from these source points
    mmp.propagate(sourcePoints);

    // Get the geodesic path from a query point to the nearest source
    SurfacePoint queryPoint2 = mesh->vertex(100);
    std::vector<SurfacePoint> path = mmp.traceBack(queryPoint2);
    Vector3 xyzCoordinates;
    std::vector<Vector3> positions;
    std::vector<std::array<size_t, 2>> edgeInds;
    double tEdgevalue;
    for (auto point:path){
        if (point.type==SurfacePointType::Vertex){
            xyzCoordinates=geometry->vertexPositions[point.vertex];
        }
        if (point.type==SurfacePointType::Edge){
            tEdgevalue=point.tEdge;
            Vertex v1=point.edge.halfedge().tailVertex();
            Vertex v2=point.edge.halfedge().tipVertex();
            xyzCoordinates=geometry->vertexPositions[v1]*tEdgevalue+
                    geometry->vertexPositions[v2]*(1-tEdgevalue);
        }
        positions.push_back(xyzCoordinates);
    }
    for(size_t i=1;i<positions.size();++i){
        edgeInds.push_back({i,i-1});
        std::cout<<positions[i]<<std::endl;
        std::cout<<edgeInds[i-1].begin()<<" "<<edgeInds[i-1].end()<<std::endl;
    }
}