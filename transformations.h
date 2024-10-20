#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include "polyhedron.h"

void rotate_point(Vertex *p, double angle, const Vertex &axis);
void translate_point(Vertex *p, double dx, double dy, double dz);
void scale_point(Vertex *p, double sx, double sy, double sz);
void reflect_point(Vertex *p, char axis);
void transform_polyhedron(Polyhedron &poly);

#endif // TRANSFORMATIONS_H