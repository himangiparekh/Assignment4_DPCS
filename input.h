#ifndef INPUT_H
#define INPUT_H

#include "polyhedron.h"

void showGuidelines();
void getInput(Polyhedron& poly);
bool validateInput(const Polyhedron& poly);
void printPolyhedron(const Polyhedron& poly);

#endif // INPUT_H