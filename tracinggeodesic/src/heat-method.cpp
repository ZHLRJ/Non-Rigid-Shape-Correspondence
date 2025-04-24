// Implement member functions HeatMethod class.
#include "heat-method.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;
    // Note: core/geometry.cpp has meanEdgeLength() function
    double meanDistance=geometry->meanEdgeLength();
    this->A = geometry->laplaceMatrix();
    this->F = geometry->massMatrix() + meanDistance*meanDistance*A ;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
//    Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
//    Vector<double> usolve = llt.solve(u);

FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {
    FaceData<Vector3> X=FaceData<Vector3>(*mesh);
    Vector<Vector3> traiangleedges(3);
    Vector3 verticeweights{0,0,0};
    for (Face f:mesh->faces()){
        Vector3 gradu{0,0,0};
        Vector3 facenormal=geometry->faceNormal(f);
        int traianglePointsIdx=0;

        for (Halfedge he:f.adjacentHalfedges()){
            traiangleedges[traianglePointsIdx]=geometry->halfedgeVector(he);
            verticeweights[(traianglePointsIdx+1)%3]=u[he.tailVertex().getIndex()];
            traianglePointsIdx+=1;
        }
        for (int i=0;i<3;i++){
            gradu+=cross(facenormal,traiangleedges[i]*verticeweights[i]);
        }
        gradu/=(2*geometry->faceArea(f));
        X[f]=-gradu.normalize();
    }

    return X; // placeholder
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {
    Vector<double> divX=Vector<double>::Zero(mesh->nVertices());
    Vector3 edge1,edge2,Xj;


    for (Vertex v:mesh->vertices()){
        double sum=0;
        for (Halfedge he:v.outgoingHalfedges()){
            Xj=X[he.face()];
            edge1=geometry->halfedgeVector(he);
            edge2=-geometry->halfedgeVector(he.next().next());
            sum+=geometry->cotan(he) *dot(edge1,Xj)+geometry->cotan(he.next().next())*dot(edge2,Xj);
        }
        divX[v.getIndex()]=sum;

    }

    return divX; // placeholder
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
    Vector<double> u = llt.solve(delta);

    FaceData<Vector3> X=computeVectorField(u);
    Vector<double> divX=computeDivergence(X);

    llt.compute(A);
    Vector<double> phi = llt.solve(-divX);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}