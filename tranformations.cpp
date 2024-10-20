#include "transformations.h"
#include <cmath>
#include <iostream>

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

            for (auto& face : poly.faces) {
                for (auto& edge : face.edges) {
                    rotate_point(&edge.v1, angle, axis);
                    rotate_point(&edge.v2, angle, axis);
                }
            }
            break;
        }
        case 2: {
            double dx, dy, dz;
            std::cout << "Enter translation values (dx dy dz): ";
            std::cin >> dx >> dy >> dz;

            for (auto& face : poly.faces) {
                for (auto& edge : face.edges) {
                    translate_point(&edge.v1, dx, dy, dz);
                    translate_point(&edge.v2, dx, dy, dz);
                }
            }
            break;
        }
        case 3: {
            double sx, sy, sz;
            std::cout << "Enter scaling factors (sx sy sz): ";
            std::cin >> sx >> sy >> sz;

            for (auto& face : poly.faces) {
                for (auto& edge : face.edges) {
                    scale_point(&edge.v1, sx, sy, sz);
                    scale_point(&edge.v2, sx, sy, sz);
                }
            }
            break;
        }
        case 4: {
            char axis;
            std::cout << "Enter axis of reflection (x, y, or z): ";
            std::cin >> axis;

            for (auto& face : poly.faces) {
                for (auto& edge : face.edges) {
                    reflect_point(&edge.v1, axis);
                    reflect_point(&edge.v2, axis);
                }
            }
            break;
        }
        default:
            std::cout << "Invalid choice\n";
    }
}