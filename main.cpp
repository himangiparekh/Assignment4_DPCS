// g++ -std=c++11 -o main main.cpp polyhedron.cpp -I /opt/homebrew/include/eigen3 -I/opt/homebrew/Cellar/sdl2/2.30.8/include -L/opt/homebrew/Cellar/sdl2/2.30.8/lib -lSDL2

#include "polyhedron.h"
#include <iostream>


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

        // Ask the user if they want to replace original input
        cout << "Do you want the original input to be replaced by the new output? (Y/N): ";
        cin.ignore();  // To consume the newline character
        cin >> replaceInput;
        if (replaceInput == 'Y' || replaceInput == 'y') {
            cout << "Original input replaced by new output.\n";
            // Implement logic to replace input with new data
        }

        // Ask if the user wants to export the results
        cout << "Do you want to export the output to a file? (Y/N): ";
        cin.ignore();  // To consume the newline character
        cin >> replaceInput;
        if (replaceInput == 'Y' || replaceInput == 'y') {
            exportOutput(poly);
        }
    }

    return 0; // Ensure main returns an integer
}