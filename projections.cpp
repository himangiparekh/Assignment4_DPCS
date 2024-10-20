#include "projections.h"
#include <SDL2/SDL.h>
#include <iostream>
#include <Eigen/Dense>

const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;

using namespace Eigen;

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
                    std::cout << "Invalid axis. Choose 'x', 'y', or 'z'.\n";
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