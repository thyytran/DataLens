#pragma once

// CFDWaveBox.h
#ifndef CFD_WAVE_BOX_H
#define CFD_WAVE_BOX_H

#include <iostream>
#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

class CFDWaveBox {
public:
    // Constructor
    CFDWaveBox(int width, int height, int depth, float gridSize);

    // Destructor
    ~CFDWaveBox();

    // Initialize OpenGL resources (VAO, VBO, shader, etc.)
    bool initializeGL();

    // Update the fluid simulation (CFD)
    void update(float deltaTime);

    // Render the fluid to the screen using OpenGL
    void render();

    // Getters for dimensions
    int getWidth() const { return width_; }
    int getHeight() const { return height_; }
    int getDepth() const { return depth_; }

private:
    // Dimensions of the simulation box (number of grid cells)
    int width_;
    int height_;
    int depth_;

    // Size of each grid cell
    float gridSize_;

    // Fluid properties (density, viscosity, etc.)
    float density_;
    float viscosity_;

    // 3D grid to store fluid data (pressure, velocity, etc.)
    std::vector<float> pressure_;
    std::vector<float> velocityX_;
    std::vector<float> velocityY_;
    std::vector<float> velocityZ_;

    // OpenGL related members
    GLuint vertexArrayObject_;
    GLuint vertexBufferObject_;
    GLuint shaderProgram_;

    // Vertex data for rendering the box (e.g., vertices of the grid)
    std::vector<float> vertices_;

    // Function to create and compile shaders
    GLuint createShader(const char* vertexShaderSource, const char* fragmentShaderSource);

    // Function to generate vertex data for the simulation box
    void generateVertices();

    // Function to apply boundary conditions
    void applyBoundaryConditions();

    // Function to calculate pressure
    void calculatePressure();

    // Function to calculate velocity
    void calculateVelocity(float deltaTime);
};

#endif // CFD_WAVE_BOX_H