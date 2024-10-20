#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <vector>
#include <Eigen/Dense>

struct Vertex {
    double x, y, z;
    bool operator<(const Vertex& other) const;
};

struct Edge {
    Vertex v1, v2;
    double length;
};

struct Face {
    std::vector<Edge> edges;
    int is_internal;
    int num_edges;
};

struct Polyhedron {
    std::vector<Face> faces;
    std::vector<Vertex> vertices;
    int num_faces;
};

struct Triangle {
    Vertex v1, v2, v3;
};

struct Tetrahedron {
    Vertex v1, v2, v3, v4;
};

// Function declarations
void showGuidelines();
void getInput(Polyhedron& poly);
void performTask(int task, Polyhedron& poly);
float calculateSurfaceArea(const Polyhedron& poly);
float calculateVolume(const Polyhedron& poly);
Vertex calculateCenterOfMass(const Polyhedron& poly);
void translatePolyhedron(Polyhedron& poly);
void rotatePolyhedron(Polyhedron& poly);
void exportOutput(Polyhedron& poly);
void orthographicProjection(Polyhedron& poly, char axis);
float calculateMomentOfInertia(const Polyhedron& poly, float mass, float inertia[3][3], const Vertex& axis);
void orthographicProjectionCustomPlane(Polyhedron& poly, float A, float B, float C, float D);
void printPolyhedron(const Polyhedron& poly);
bool validateInput(const Polyhedron& poly);
void transform_polyhedron(Polyhedron& poly);

// Helper function declarations
Vertex vectorSubtract(const Vertex& a, const Vertex& b);
float vectorDot(const Vertex& a, const Vertex& b);
Vertex vectorCross(const Vertex& a, const Vertex& b);
float vectorMagnitude(const Vertex& v);
float calculateTetrahedronVolume(const Vertex& v0, const Vertex& v1, const Vertex& v2);
Vertex calculateTetrahedronCentroid(const Vertex& v0, const Vertex& v1, const Vertex& v2);
void rotate_point(Vertex* p, double angle, const Vertex& axis);
void translate_point(Vertex* p, double dx, double dy, double dz);
void scale_point(Vertex* p, double sx, double sy, double sz);
void reflect_point(Vertex* p, char axis);

#endif // POLYHEDRON_H