// g++ -std=c++11 -o main2 tryfull6.cpp -I /opt/homebrew/include/eigen3 -I/opt/homebrew/Cellar/sdl2/2.30.8/include -L/opt/homebrew/Cellar/sdl2/2.30.8/lib -lSDL2

// FINAL: RECONSTRUCTION + PROJECTION + SLICING + GEOMETRIC + TRANSFORMATION

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <SDL2/SDL.h>
#include <Eigen/Dense>
#include <map>
#include <algorithm> 

using namespace Eigen;
using namespace std;

// Constants for SDL Window dimensions
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;

struct Vertex {
    
    double x, y, z;

    // Define the '<' operator to compare Vertex objects
    bool operator<(const Vertex& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};

struct Edge {
    Vertex v1, v2;
    double length;
};

struct Face {
    vector<Edge> edges;
    int is_internal;
    int num_edges; 
};

struct Polyhedron {
    vector<Face> faces;
    vector<Vertex> vertices;
    int num_faces;
};

struct Triangle {
    Vertex v1, v2, v3;
};

struct Tetrahedron {
    Vertex v1, v2, v3, v4;
};

// Function prototypes
void showGuidelines();
void getInput(Polyhedron& poly);
float calculateSurfaceArea(const Polyhedron& poly);
float calculateVolume(const Polyhedron& poly);
Vertex calculateCenterOfMass(const Polyhedron& poly);
void translatePolyhedron(Polyhedron& poly);
void rotatePolyhedron(Polyhedron& poly);
void exportOutput(Polyhedron& poly);
void orthographicProjection(Polyhedron& poly, char axis); // done
float calculateMomentOfInertia(const Polyhedron& poly, float mass, float inertia[3][3], const Vertex& axis);
void orthographicProjectionCustomPlane(Polyhedron& poly, float A, float B, float C, float D); // done
void printPolyhedron(const Polyhedron& poly); // done
bool validateInput(const Polyhedron& poly); // done

// Helper functions for vector math
Vertex vectorSubtract(const Vertex& a, const Vertex& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

float vectorDot(const Vertex& a, const Vertex& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vertex vectorCross(const Vertex& a, const Vertex& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

float vectorMagnitude(const Vertex& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Function to calculate the signed volume of a tetrahedron
float calculateTetrahedronVolume(const Vertex& v0, const Vertex& v1, const Vertex& v2) {
    return vectorDot(v0, vectorCross(v1, v2)) / 6.0f;
}

// Function to calculate the centroid of a tetrahedron
Vertex calculateTetrahedronCentroid(const Vertex& v0, const Vertex& v1, const Vertex& v2) {
    return {
        (v0.x + v1.x + v2.x) / 4.0f,
        (v0.y + v1.y + v2.y) / 4.0f,
        (v0.z + v1.z + v2.z) / 4.0f
    };
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

float calculateAxisMomentOfInertia(const float inertia[3][3], const Vertex& axis) {
    float Iaxis = 0.0f;
    Iaxis += inertia[0][0] * axis.x * axis.x;
    Iaxis += inertia[1][1] * axis.y * axis.y;
    Iaxis += inertia[2][2] * axis.z * axis.z;
    Iaxis += 2 * inertia[0][1] * axis.x * axis.y;
    Iaxis += 2 * inertia[0][2] * axis.x * axis.z;
    Iaxis += 2 * inertia[1][2] * axis.y * axis.z;
    return Iaxis;
}

float calculateMomentOfInertia(const Polyhedron& poly, float mass, float inertia[3][3], const Vertex& axis) {
    Vertex com = calculateCenterOfMass(poly); // Corrected function call
    float volume = calculateVolume(poly);
    float density = mass / volume;

    // Normalize the axis vector
    float axisLength = std::sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
    Vertex normalizedAxis = {axis.x / axisLength, axis.y / axisLength, axis.z / axisLength};

    // Initialize inertia matrix to zero
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inertia[i][j] = 0.0f;
        }
    }

    // Calculate the inertia tensor
    for (const auto& face : poly.faces) {
        if (face.edges.size() < 3) continue;

        Vertex v0 = vectorSubtract(face.edges[0].v1, com);

        for (size_t e = 1; e < face.edges.size() - 1; e++) {
            Vertex v1 = vectorSubtract(face.edges[e].v1, com);
            Vertex v2 = vectorSubtract(face.edges[e + 1].v1, com);

            Vertex crossProd = vectorCross(v1, v2);
            float detJ = vectorDot(v0, crossProd);

            float subexpr0 = vectorDot(v0, v0) + vectorDot(v1, v1) + vectorDot(v2, v2);
            float subexpr1 = v0.x * v1.x + v0.y * v1.y + v0.z * v1.z + 
                             v0.x * v2.x + v0.y * v2.y + v0.z * v2.z + 
                             v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
            float subexpr2 = v0.x + v1.x + v2.x;
            float subexpr3 = v0.y + v1.y + v2.y;

            inertia[0][0] += detJ * (subexpr0 + subexpr1);
            inertia[1][1] += detJ * (subexpr0 + subexpr1);
            inertia[2][2] += detJ * (subexpr0 + subexpr1);
            inertia[0][1] -= detJ * subexpr2 * subexpr3;
            inertia[0][2] -= detJ * subexpr2 * (v0.z + v1.z + v2.z);
            inertia[1][2] -= detJ * subexpr3 * (v0.z + v1.z + v2.z);
        }
    }

    // Scale the inertia tensor by density
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inertia[i][j] *= density / 60.0f; // Assuming a constant scaling factor
            if (i != j) {
                inertia[j][i] = inertia[i][j]; // Make the inertia matrix symmetric
            }
        }
    }

    // Return the moment of inertia about the specified axis
    return calculateAxisMomentOfInertia(inertia, normalizedAxis);
}


// Surface Area calculation
float calculateSurfaceArea(const Polyhedron& poly) {
    float totalArea = 0.0f;

    for (const auto& face : poly.faces) {
        if (face.edges.size() < 3 || face.is_internal) continue;

        Vertex v0 = face.edges[0].v1;

        for (size_t j = 1; j < face.edges.size() - 1; j++) {
            Vertex v1 = face.edges[j].v1;
            Vertex v2 = face.edges[j + 1].v1;

            Vertex edge1 = vectorSubtract(v1, v0);
            Vertex edge2 = vectorSubtract(v2, v0);

            Vertex crossProd = vectorCross(edge1, edge2);
            float triangleArea = vectorMagnitude(crossProd) / 2.0f;

            totalArea += triangleArea;
        }
    }

    return totalArea;
}

// Volume calculation
float calculateVolume(const Polyhedron& poly) {
    float totalVolume = 0.0f;

    for (const auto& face : poly.faces) {
        if (face.is_internal) continue;

        Vertex v0 = face.edges[0].v1;

        for (size_t j = 1; j < face.edges.size() - 1; j++) {
            Vertex v1 = face.edges[j].v1;
            Vertex v2 = face.edges[j + 1].v1;

            totalVolume += calculateTetrahedronVolume(v0, v1, v2);
        }
    }

    return totalVolume;
}

Vertex projectVertexOntoPlane(Vertex v, float A, float B, float C, float D) {
    Vertex projected;
    float dotProduct = v.x * A + v.y * B + v.z * C;  // P · n
    float normalLenSquared = A * A + B * B + C * C;  // n · n
    float scale = (dotProduct + D) / normalLenSquared;

    // Projected point: P' = P - scale * n
    projected.x = v.x - scale * A;
    projected.y = v.y - scale * B;
    projected.z = v.z - scale * C;

    return projected;
}

float dotProduct(float x1, float y1, float z1, float x2, float y2, float z2) {
    return x1 * x2 + y1 * y2 + z1 * z2;
}

// Function to normalize a 3D vector
void normalizeVector(float* x, float* y, float* z) {
    float length = sqrt(*x * *x + *y * *y + *z * *z);
    *x /= length;
    *y /= length;
    *z /= length;
}

// Function to print the reconstructed 3D polyhedron
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

// Helper function to check collinearity of three points
bool checkCollinearity(const Vertex& p1, const Vertex& p2, const Vertex& p3) {
    Vector3d v1(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    Vector3d v2(p3.x - p2.x, p3.y - p2.y, p3.z - p2.z);
    Vector3d crossProduct = v1.cross(v2);
    return crossProduct.norm() < 1e-6;  // If the cross product is near zero, points are collinear
}

// Helper function to check if a set of points lie on the same plane (planarity)
bool checkPlanarity(const vector<Vertex>& face) {
    if (face.size() < 3) return true; // Less than 3 points always planar
    Vector3d normal(0, 0, 0);
    for (size_t i = 1; i < face.size() - 1; ++i) {
        Vector3d v1(face[i].x - face[0].x, face[i].y - face[0].y, face[i].z - face[0].z);
        Vector3d v2(face[i + 1].x - face[0].x, face[i + 1].y - face[0].y, face[i + 1].z - face[0].z);
        normal += v1.cross(v2);
    }
    return normal.norm() > 1e-6;  // If the normal is near zero, the points are not planar
}

// Helper function to compute distance between two points
double computeDistance(const Vertex& p1, const Vertex& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

// Function to validate the input projections and 3D reconstruction
bool validateInput(const Polyhedron& poly) {
    // 1. Edge length consistency
    for (size_t i = 0; i < poly.faces.size(); i++) {
        for (size_t j = 0; j < poly.faces[i].edges.size(); j++) {
            const Edge& edge = poly.faces[i].edges[j];
            double dist3D = computeDistance(edge.v1, edge.v2);
            if (fabs(dist3D - edge.length) > 1e-6) {
                printf("Edge length inconsistency detected for edge %zu in face %zu\n", j + 1, i + 1);
                return false;
            }
        }
    }

    // 2. Collinearity check and 3. Planarity check for faces
    for (size_t i = 0; i < poly.faces.size(); i++) {
        vector<Vertex> faceVertices;
        for (size_t j = 0; j < poly.faces[i].edges.size(); j++) {
            faceVertices.push_back(poly.faces[i].edges[j].v1);
        }
        
        if (faceVertices.size() >= 3) {
            for (size_t k = 0; k < faceVertices.size() - 2; ++k) {
                if (checkCollinearity(faceVertices[k], faceVertices[k + 1], faceVertices[k + 2])) {
                    printf("Collinearity detected for points in face %zu\n", i + 1);
                    return false;
                }
            }
        }
        
        if (!checkPlanarity(faceVertices)) {
            printf("Non-planar face detected for face %zu\n", i + 1);
            return false;
        }
    }

    // 4. Check if the polyhedron is closed
    map<pair<Vertex, Vertex>, int> edgeCount;
    
    for (const auto& face : poly.faces) {
        for (const auto& edge : face.edges) {
            // Use an unordered pair of vertices to represent an edge, ensuring consistent ordering
            pair<Vertex, Vertex> edgeKey = minmax(edge.v1, edge.v2);
            edgeCount[edgeKey]++;
        }
    }

// Each edge should appear exactly twice (once for each adjacent face)
    for (const auto& entry : edgeCount) {
        if (entry.second != 2) {
            // Print the coordinates of the two vertices
            printf("Edge between vertices (%.2f, %.2f, %.2f) and (%.2f, %.2f, %.2f) is not shared by exactly two faces\n",
                entry.first.first.x, entry.first.first.y, entry.first.first.z,
                entry.first.second.x, entry.first.second.y, entry.first.second.z);
            return false;
        }
    }

    return true;
}

// Show user guidelines at the start
void showGuidelines() {
    cout << "Here are the specifications and user guidelines.\n";
    cout << "You will be prompted to provide mandatory inputs.\n";
    cout << "You can choose tasks such as geometric calculations, projections, or transformations.\n";
    cout << "At the end of each task, you can decide if the original input should be replaced.\n";
    cout << "Use the command-line interface to interact with the program.\n";
    cout << "At any point, you can choose to export the results to a file.\n";
}

// Generalized Center of Mass calculation using tetrahedral decomposition
Vertex calculateCenterOfMass(const Polyhedron& poly) {
    float totalVolume = 0.0f;
    Vertex tempCenterOfMass = {0, 0, 0};

    for (const auto& face : poly.faces) {
        if (face.edges.size() < 3 || face.is_internal) continue;

        Vertex v0 = face.edges[0].v1;

        for (size_t j = 1; j < face.edges.size() - 1; j++) {
            Vertex v1 = face.edges[j].v1;
            Vertex v2 = face.edges[j + 1].v1;

            float tetrahedronVolume = calculateTetrahedronVolume(v0, v1, v2);
            Vertex tetrahedronCentroid = calculateTetrahedronCentroid(v0, v1, v2);

            tempCenterOfMass.x += tetrahedronCentroid.x * tetrahedronVolume;
            tempCenterOfMass.y += tetrahedronCentroid.y * tetrahedronVolume;
            tempCenterOfMass.z += tetrahedronCentroid.z * tetrahedronVolume;

            totalVolume += tetrahedronVolume;
        }
    }

    if (totalVolume > 0) {
        tempCenterOfMass.x /= totalVolume;
        tempCenterOfMass.y /= totalVolume;
        tempCenterOfMass.z /= totalVolume;
    }

    return tempCenterOfMass;
}

void orthographicProjection(Polyhedron& poly, char axis) {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow("Orthographic Projection", 
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    // Loop through faces and project each edge
    for (const auto& face : poly.faces) {
        for (const auto& edge : face.edges) {
            int x1, y1, x2, y2;

            switch (axis) {
                case 'x':  // Project onto the YZ plane
                    x1 = static_cast<int>(edge.v1.y * 100 + WINDOW_WIDTH / 2);
                    y1 = static_cast<int>(edge.v1.z * 100 + WINDOW_HEIGHT / 2);
                    x2 = static_cast<int>(edge.v2.y * 100 + WINDOW_WIDTH / 2);
                    y2 = static_cast<int>(edge.v2.z * 100 + WINDOW_HEIGHT / 2);
                    break;
                case 'y':  // Project onto the XZ plane
                    x1 = static_cast<int>(edge.v1.x * 100 + WINDOW_WIDTH / 2);
                    y1 = static_cast<int>(edge.v1.z * 100 + WINDOW_HEIGHT / 2);
                    x2 = static_cast<int>(edge.v2.x * 100 + WINDOW_WIDTH / 2);
                    y2 = static_cast<int>(edge.v2.z * 100 + WINDOW_HEIGHT / 2);
                    break;
                case 'z':  // Project onto the XY plane
                    x1 = static_cast<int>(edge.v1.x * 100 + WINDOW_WIDTH / 2);
                    y1 = static_cast<int>(edge.v1.y * 100 + WINDOW_HEIGHT / 2);
                    x2 = static_cast<int>(edge.v2.x * 100 + WINDOW_WIDTH / 2);
                    y2 = static_cast<int>(edge.v2.y * 100 + WINDOW_HEIGHT / 2);
                    break;
                default:
                    cout << "Invalid axis. Choose 'x', 'y', or 'z'.\n";
                    return;
            }

            SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
        }
    }

    SDL_RenderPresent(renderer);  // Present the rendered image to the screen

    // Wait for a quit event to close the window
    SDL_Event event;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }
        }
    }

    // Cleanup
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void orthographicProjectionCustomPlane(Polyhedron& poly, float A, float B, float C, float D) {
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window* window = SDL_CreateWindow("Orthographic Projection on Custom Plane", 
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_SHOWN);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    // Loop through faces and project each edge
    for (const auto& face : poly.faces) {
        for (const auto& edge : face.edges) {
            // Determine the color based on the edge's orientation
            Vector3d v1(edge.v1.x, edge.v1.y, edge.v1.z);
            Vector3d v2(edge.v2.x, edge.v2.y, edge.v2.z);
            Vector3d edgeDirection = v2 - v1;

            // Normal vector of the plane
            Vector3d planeNormal(A, B, C);
            double d = D;  // Plane equation constant

            // Calculate the dot product to determine facing direction
            if (edgeDirection.dot(planeNormal) > 0) {
                SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);  // Front-facing edges are black
            } else {
                SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);  // Back-facing edges are red
            }

            // Project the edge onto the plane
            // Simple orthographic projection ignoring the normal direction
            int x1 = static_cast<int>(edge.v1.x * 100 + WINDOW_WIDTH / 2);
            int y1 = static_cast<int>(edge.v1.y * 100 + WINDOW_HEIGHT / 2);
            int x2 = static_cast<int>(edge.v2.x * 100 + WINDOW_WIDTH / 2);
            int y2 = static_cast<int>(edge.v2.y * 100 + WINDOW_HEIGHT / 2);

            SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
        }
    }

    SDL_RenderPresent(renderer);  // Present the rendered image to the screen

    // Wait for a quit event to close the window
    SDL_Event event;
    bool quit = false;
    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }
        }
    }

    // Cleanup
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void rotate_point(Vertex *p, double angle, const Vertex &axis) {
    double rad = angle * M_PI / 180.0;
    double c = cos(rad);
    double s = sin(rad);
    double t = 1.0 - c;

    // Normalize the axis vector
    double magnitude = sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
    Vertex normAxis = {axis.x / magnitude, axis.y / magnitude, axis.z / magnitude};

    // Rotation matrix elements
    double m11 = t * normAxis.x * normAxis.x + c;
    double m12 = t * normAxis.x * normAxis.y - s * normAxis.z;
    double m13 = t * normAxis.x * normAxis.z + s * normAxis.y;
    double m21 = t * normAxis.x * normAxis.y + s * normAxis.z;
    double m22 = t * normAxis.y * normAxis.y + c;
    double m23 = t * normAxis.y * normAxis.z - s * normAxis.x;
    double m31 = t * normAxis.x * normAxis.z - s * normAxis.y;
    double m32 = t * normAxis.y * normAxis.z + s * normAxis.x;
    double m33 = t * normAxis.z * normAxis.z + c;

    // Apply rotation
    double x = p->x, y = p->y, z = p->z;
    p->x = m11 * x + m12 * y + m13 * z;
    p->y = m21 * x + m22 * y + m23 * z;
    p->z = m31 * x + m32 * y + m33 * z;
}

void translate_point(Vertex *p, double dx, double dy, double dz) {
    p->x += dx;
    p->y += dy;
    p->z += dz;
}

void scale_point(Vertex *p, double sx, double sy, double sz) {
    p->x *= sx;
    p->y *= sy;
    p->z *= sz;
}

void reflect_point(Vertex *p, char axis) {
    switch (axis) {
        case 'x': p->x = -p->x; break;
        case 'y': p->y = -p->y; break;
        case 'z': p->z = -p->z; break;
    }
}

void transform_polyhedron(Polyhedron &poly) {
    int choice;
    std::cout << "Choose transformation:\n";
    std::cout << "1. Rotate\n2. Translate\n3. Scale\n4. Reflect\n";
    std::cin >> choice;

    switch (choice) {
        case 1: {
            double angle;
            Vertex axis;
            std::cout << "Enter rotation angle (in degrees): ";
            std::cin >> angle;
            std::cout << "Enter rotation axis (x y z): ";
            std::cin >> axis.x >> axis.y >> axis.z;

            std::cout << "\nRotated Polyhedron:\n";
            for (int i = 0; i < poly.num_faces; i++) {
                for (int j = 0; j < poly.faces[i].num_edges; j++) {
                    rotate_point(&poly.faces[i].edges[j].v1, angle, axis);
                    rotate_point(&poly.faces[i].edges[j].v2, angle, axis);
                }
            }
            printPolyhedron(poly);
            break;
        }
        case 2: {
            double dx, dy, dz;
            std::cout << "Enter translation values (dx dy dz): ";
            std::cin >> dx >> dy >> dz;

            std::cout << "\nTranslated Polyhedron:\n";
            for (int i = 0; i < poly.num_faces; i++) {
                for (int j = 0; j < poly.faces[i].num_edges; j++) {
                    translate_point(&poly.faces[i].edges[j].v1, dx, dy, dz);
                    translate_point(&poly.faces[i].edges[j].v2, dx, dy, dz);
                }
            }
            printPolyhedron(poly);
            break;
        }
        case 3: {
            double sx, sy, sz;
            std::cout << "Enter scaling factors (sx sy sz): ";
            std::cin >> sx >> sy >> sz;

            std::cout << "\nScaled Polyhedron:\n";
            for (int i = 0; i < poly.num_faces; i++) {
                for (int j = 0; j < poly.faces[i].num_edges; j++) {
                    scale_point(&poly.faces[i].edges[j].v1, sx, sy, sz);
                    scale_point(&poly.faces[i].edges[j].v2, sx, sy, sz);
                }
            }
            printPolyhedron(poly);
            break;
        }
        case 4: {
            char axis;
            std::cout << "Enter axis of reflection (x, y, or z): ";
            std::cin >> axis;

            std::cout << "\nReflected Polyhedron:\n";
            for (int i = 0; i < poly.num_faces; i++) {
                for (int j = 0; j < poly.faces[i].num_edges; j++) {
                    reflect_point(&poly.faces[i].edges[j].v1, axis);
                    reflect_point(&poly.faces[i].edges[j].v2, axis);
                }
            }
            printPolyhedron(poly);
            break;
        }
        default:
            std::cout << "Invalid choice\n";
    }
}


int main() {
    Polyhedron poly;
    int task;
    char replaceInput;
    float mass;
    Vertex axis;
    float inertiaMatrix[3][3]; // Changed name to avoid conflict
    float momentOfInertia; // Variable to hold the calculated moment of inertia

    showGuidelines();  // Show specifications and guidelines to the user

    // Get the mandatory input from the user
    getInput(poly);

    // Print the reconstructed polyhedron
    printPolyhedron(poly);

    // Validate the input
    bool isValid = validateInput(poly);
    
    if (isValid) {
        cout << "The input and reconstruction are valid.\n";
    } else {
        cout << "The input or reconstruction is invalid.\n";
        return 1;  // Exit if the input is invalid
    }

    // Main loop for performing tasks
    while (true) {
        cout << "\nSelect a task to perform:\n";
        cout << "1. Calculate Surface Area\n";
        cout << "2. Calculate Volume\n";
        cout << "3. Calculate Center of Mass\n";
        cout << "4. Transform Polyhedron\n";
        cout << "6. Orthographic Projection\n";
        cout << "7. Orthographic Projection onto Custom Plane\n";
        cout << "8. Calculate Moment of Inertia\n";
        cout << "9. Exit\n";
        cout << "Enter your choice: ";
        cin >> task;

        if (task == 9) {
            cout << "Exiting the program.\n";
            break;
        }

        if (task == 4) {  // Transform Polyhedron
            transform_polyhedron(poly);
        }

        if (task == 1) {  // Calculate Surface Area
            float surfaceArea = calculateSurfaceArea(poly);
            cout << "The Surface Area of the polyhedron is: " << surfaceArea << endl;
        }

        if (task == 2) {  // Calculate Volume
            float volume = calculateVolume(poly);
            cout << "The Volume of the polyhedron is: " << volume << endl;
        }
        
        if (task == 3) {  // Calculate Centre of Mass
            Vertex centre = calculateCenterOfMass(poly);
            cout << "The Centre of Mass of the polyhedron is: " << centre.x << ", " << centre.y << ", " << centre.z << endl;
        }

        if (task == 8) {  // Calculate Moment of Inertia
            std::cout << "Enter the mass of the polyhedron: ";
            std::cin >> mass;
            if (mass <= 0) {
                std::cerr << "Mass must be a positive value." << std::endl;
                return 1; // Exit the program with an error code
            }

            std::cout << "Enter the axis for moment of inertia calculation (x y z): Put in co-efficients for the axis: ";
            std::cin >> axis.x >> axis.y >> axis.z;

            momentOfInertia = calculateMomentOfInertia(poly, mass, inertiaMatrix, axis); // Call function
            cout << "The Moment of Inertia of the polyhedron about the specified axis is: " << momentOfInertia << endl; // Correct output message
        }

        if (task == 6) {
            char axis;
            cout << "Select axis for orthographic projection (x, y, z): ";
            cin.ignore();  // Consume the newline character
            cin >> axis;
            orthographicProjection(poly, axis);
        }

        if (task == 7) {
            float A, B, C, D;
            std::cout << "Enter the coefficients of the plane equation (Ax + By + Cz = D):" << std::endl;
            std::cout << "A (normal x-component): ";
            std::cin >> A;
            std::cout << "B (normal y-component): ";
            std::cin >> B;
            std::cout << "C (normal z-component): ";
            std::cin >> C;
            std::cout << "D (distance from origin): ";
            std::cin >> D;
            orthographicProjectionCustomPlane(poly, A, B, C, D);
        }
    }

    return 0; // Ensure main returns an integer
}