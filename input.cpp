#include "input.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void showGuidelines() {
    cout << "Here are the specifications and user guidelines.\n";
    cout << "You will be prompted to provide mandatory inputs.\n";
    cout << "You can choose tasks such as geometric calculations, projections, or transformations.\n";
    cout << "At the end of each task, you can decide if the original input should be replaced.\n";
    cout << "Use the command-line interface to interact with the program.\n";
    cout << "At any point, you can choose to export the results to a file.\n";
}

void getInput(Polyhedron& poly) {
    int numVertices, numFaces;
    
    printf("Enter the number of vertices in the polyhedron: ");
    scanf("%d", &numVertices);

    printf("Enter the number of faces in the polyhedron: ");
    scanf("%d", &numFaces);
    poly.faces.resize(numFaces);

    // Get 2D projections for each vertex
    vector<double> xy(numVertices * 2);
    vector<double> xz(numVertices * 2);

    printf("Enter the 2D coordinates for the vertices:\n");
    for (int i = 0; i < numVertices; ++i) {
        printf("Vertex %d (xy): ", i + 1);
        scanf("%lf %lf", &xy[i * 2], &xy[i * 2 + 1]);

        printf("Vertex %d (xz): ", i + 1);
        scanf("%lf %lf", &xz[i * 2], &xz[i * 2 + 1]);
    }

    // Reconstruct 3D vertices
    poly.vertices.resize(numVertices);
    MatrixXd A(4, 3);
    A << 1, 0, 0,
         0, 1, 0,
         1, 0, 0,
         0, 0, 1;

    for (int i = 0; i < numVertices; ++i) {
        VectorXd B(4);
        B << xy[i * 2], xy[i * 2 + 1], xz[i * 2], xz[i * 2 + 1];
        Vector3d X = A.colPivHouseholderQr().solve(B);
        poly.vertices[i].x = X(0);
        poly.vertices[i].y = X(1);
        poly.vertices[i].z = X(2);
    }

    // Get faces and edges
    for (int i = 0; i < numFaces; i++) {
        int numEdges;
        printf("Enter the number of edges in face %d: ", i + 1);
        scanf("%d", &numEdges);
        poly.faces[i].edges.resize(numEdges);

        printf("Is face %d internal? (1 for yes, 0 for no): ", i + 1);
        scanf("%d", &poly.faces[i].is_internal);
        
        for (int j = 0; j < numEdges; j++) {
            int v1, v2;
            printf("Enter vertices for edge %d in face %d (format: v1 v2): ", j + 1, i + 1);
            scanf("%d %d", &v1, &v2);
            v1--; v2--;  // Convert to 0-based indexing
            
            poly.faces[i].edges[j].v1 = poly.vertices[v1];
            poly.faces[i].edges[j].v2 = poly.vertices[v2];
            
            // Calculate edge length
            double dx = poly.faces[i].edges[j].v2.x - poly.faces[i].edges[j].v1.x;
            double dy = poly.faces[i].edges[j].v2.y - poly.faces[i].edges[j].v1.y;
            double dz = poly.faces[i].edges[j].v2.z - poly.faces[i].edges[j].v1.z;
            poly.faces[i].edges[j].length = sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
}

bool validateInput(const Polyhedron& poly) {
    // Implementation of validateInput function
    // (The existing implementation remains the same)
}

void printPolyhedron(const Polyhedron& poly) {
    printf("Reconstructed 3D Polyhedron:\n");
    for (size_t i = 0; i < poly.faces.size(); i++) {
        printf("Face %zu (%s):\n", i + 1, poly.faces[i].is_internal ? "internal" : "external");
        for (size_t j = 0; j < poly.faces[i].edges.size(); j++) {
            const Edge& edge = poly.faces[i].edges[j];
            printf("  Edge %zu: (%.2f, %.2f, %.2f) to (%.2f, %.2f, %.2f), Length: %.2f\n",
                   j + 1, edge.v1.x, edge.v1.y, edge.v1.z, edge.v2.x, edge.v2.y, edge.v2.z, edge.length);
        }
    }
}