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

//VertexData<int>* dist;
VertexData<int> dist;
std::vector<int> shuffleVector(int N,bool fixedSeed=true){
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
int main(){
    std::vector<std::array<double, 3>> colorbank{ {1,0,0},
                                                  {0,1,0},
                                                  {0,0,1}};

    auto start = std::chrono::steady_clock::now();
    for (auto item:colorbank){
        std::cout<<item[0]<<std::endl;
    }
    auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout << "It took " << elapsed_seconds.count() << "s" << std::endl;
//    for (std::vector<float>::iterator item=v.begin();item<v.end();++item){
//        std::cout<<*item<<std::endl;
////        Address of iterator can be printed by:
//        std::cout << &item <<std:: endl;
////        Address of iterated elements can be printed by:
//        std::cout << &(*item) <<std:: endl;
//
//    }




//    std::vector<int> v (5);
//    std::iota(v.begin(),v.end(),4);
//    // Show the picked vertex.
//    //    std::vector<Vector3> vertPos;
//    //    std::vector<std::array<size_t, 2>> vertInd;
//    for (std::vector<int>::iterator it = v.begin(); it != v.end(); ++it) {
//        std::cout<<*it<<std::endl;
//
////        vertPos.push_back(geometry->inputVertexPositions[*it]);
//    }


//    std::random_device rd;
//    std::mt19937 g(0);
//
//    // Use of shuffle
////    std::shuffle(v.begin(), v.end(), default_random_engine(seed));
//    std::shuffle(v.begin(), v.end(), g);
//
//    std::copy(v.begin(), v.end(), std::ostream_iterator<int>(std::cout, " "));
//    std::cout << "\n";
//    srand(17):
//    Iteration 1 ranvals [0] = 285719
//    Iteration 2 ranvals [1] = 507111939
//    Iteration 3 ranvals [2] = 1815247477
//    Iteration 4 ranvals [3] = 1711656657
//    Iteration 5 ranvals [4] = 122498987
//    srand(1):
//    Iteration 1 ranvals [0] = 16807
//    Iteration 2 ranvals [1] = 282475249
//    Iteration 3 ranvals [2] = 1622650073
//    Iteration 4 ranvals [3] = 984943658
//    Iteration 5 ranvals [4] = 1144108930
//    std::string filepath = "../../../../input/bunny.obj";
//    std::tie(mesh_uptr, geometry_uptr) = readManifoldSurfaceMesh(filepath);
//    geometry=geometry_uptr.release();
//    mesh=mesh_uptr.release();
//
//    VertexData<double> dist(*mesh, 0.);
//    VectorHeatMethodSolver vhtSolver(*geometry);
//
//    SurfacePoint center = findCenter(*mesh,*geometry, vhtSolver, dist, 2);


}