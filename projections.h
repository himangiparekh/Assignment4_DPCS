#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include "polyhedron.h"

void orthographicProjection(Polyhedron& poly, char axis);
void orthographicProjectionCustomPlane(Polyhedron& poly, float A, float B, float C, float D);

#endif // PROJECTIONS_H