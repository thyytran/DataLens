/*
// CFDWaveBox.cpp
#include "include/cfd/CFDWaveBox.h"
#include <cmath> // For math functions like sin, cos

// Simple vertex shader source code
const char* vertexShaderSource = R"(
    #version 330 core
    layout (location = 0) in vec3 aPos;

    void main() {
        gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);
    }
)";

// Simple fragment shader source code
const char* fragmentShaderSource = R"(
    #version 330 core
    out vec4 FragColor;

    void main() {
        FragColor = vec4(0.0, 0.5, 1.0, 1.0); // Blue color
    }
)";


CFDWaveBox::CFDWaveBox(int width, int height, int depth, float gridSize)
    : width_(width), height_(height), depth_(depth), gridSize_(gridSize),
    density_(1.0f), viscosity_(0.01f), vertexArrayObject_(0), vertexBufferObject_(0), shaderProgram_(0) {

    // Resize the fluid data vectors
    int totalCells = width_ * height_ * depth_;
    pressure_.resize(totalCells, 0.0f);
    velocityX_.resize(totalCells, 0.0f);
    velocityY_.resize(totalCells, 0.0f);
    velocityZ_.resize(totalCells, 0.0f);

    // Initialize some wave source (example: a sine wave at the top surface)
    for (int x = 0; x < width_; ++x) {
        for (int z = 0; z < depth_; ++z) {
            int index = x + z * width_ + (height_ - 1) * width_ * depth_; // Top layer
            pressure_[index] = 0.5f * sin(x * 0.2f) * cos(z * 0.2f);
        }
    }
}

CFDWaveBox::~CFDWaveBox() {
    // Cleanup OpenGL resources
    glDeleteVertexArrays(1, &vertexArrayObject_);
    glDeleteBuffers(1, &vertexBufferObject_);
    glDeleteProgram(shaderProgram_);
}

bool CFDWaveBox::initializeGL() {
    // Generate vertices for the simulation box
    generateVertices();

    // Create and compile shaders
    shaderProgram_ = createShader(vertexShaderSource, fragmentShaderSource);
    if (shaderProgram_ == 0) {
        std::cerr << "Err > Failed to create shader program.\n\n";
        return false;
    }

    // Generate vertex array object (VAO)
    glGenVertexArrays(1, &vertexArrayObject_);
    glBindVertexArray(vertexArrayObject_);

    // Generate vertex buffer object (VBO)
    glGenBuffers(1, &vertexBufferObject_);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferObject_);
    glBufferData(GL_ARRAY_BUFFER, vertices_.size() * sizeof(float), vertices_.data(), GL_STATIC_DRAW);

    // Set vertex attributes
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Unbind VAO and VBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    return true;
}

void CFDWaveBox::update(float deltaTime) {
    // 1. Apply boundary conditions
    applyBoundaryConditions();

    // 2. Calculate pressure
    calculatePressure();

    // 3. Calculate velocity
    calculateVelocity(deltaTime);
}

void CFDWaveBox::render() {
    // Use the shader program
    glUseProgram(shaderProgram_);

    // Bind the VAO
    glBindVertexArray(vertexArrayObject_);

    // Draw the simulation box (as points for now)
    glDrawArrays(GL_POINTS, 0, vertices_.size() / 3);  // Assuming 3 floats per vertex

    // Unbind the VAO
    glBindVertexArray(0);
}


GLuint CFDWaveBox::createShader(const char* vertexShaderSource, const char* fragmentShaderSource) {
    // Vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // Check vertex shader compilation errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "Err > Vertex shader compilation failed:\n" << infoLog << "\n\n";
        return 0;
    }

    // Fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // Check fragment shader compilation errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "Err > Fragment shader compilation failed:\n" << infoLog << "\n\n";
        glDeleteShader(vertexShader);
        return 0;
    }

    // Shader program
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // Check shader program linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "Err > Shader program linking failed:\n" << infoLog << "\n\n";
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        return 0;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return shaderProgram;
}


void CFDWaveBox::generateVertices() {
    // Clear existing vertices
    vertices_.clear();

    // Generate vertices for each grid cell
    for (int x = 0; x < width_; ++x) {
        for (int y = 0; y < height_; ++y) {
            for (int z = 0; z < depth_; ++z) {
                // Calculate the position of the grid cell
                float posX = x * gridSize_;
                float posY = y * gridSize_;
                float posZ = z * gridSize_;

                // Add the vertex to the vertices vector (just a point for now)
                vertices_.push_back((posX / (width_ * gridSize_)) * 2.0f - 1.0f);   // Normalize to -1 to 1 range
                vertices_.push_back((posY / (height_ * gridSize_)) * 2.0f - 1.0f);  // Normalize to -1 to 1 range
                vertices_.push_back((posZ / (depth_ * gridSize_)) * 2.0f - 1.0f);   // Normalize to -1 to 1 range
            }
        }
    }
}


void CFDWaveBox::applyBoundaryConditions() {
    // Example: Set pressure to 0 at the boundaries
    for (int x = 0; x < width_; ++x) {
        for (int z = 0; z < depth_; ++z) {
            pressure_[x + z * width_] = 0.0f; // Bottom
            pressure_[x + z * width_ + (height_ - 1) * width_ * depth_] = 0.0f; // Top
        }
    }
    for (int y = 0; y < height_; ++y) {
        for (int z = 0; z < depth_; ++z) {
            pressure_[z * width_ * depth_ + y * width_] = 0.0f; // Left
            pressure_[(width_ - 1) + z * width_ * depth_ + y * width_] = 0.0f; // Right
        }
    }
    for (int x = 0; x < width_; ++x) {
        for (int y = 0; y < height_; ++y) {
            pressure_[x + y * width_] = 0.0f; // Front
            pressure_[x + y * width_ + (depth_ - 1) * width_ * depth_] = 0.0f; // Back
        }
    }
}

void CFDWaveBox::calculatePressure() {
    // Simple pressure calculation (example: average of neighbors)
    std::vector<float> newPressure_ = pressure_; // Create a copy
    for (int x = 1; x < width_ - 1; ++x) {
        for (int y = 1; y < height_ - 1; ++y) {
            for (int z = 1; z < depth_ - 1; ++z) {
                int index = x + y * width_ * depth_ + z * width_;
                newPressure_[index] = (
                    pressure_[index - 1] +   // Left
                    pressure_[index + 1] +   // Right
                    pressure_[index - width_ * depth_] + // Below
                    pressure_[index + width_ * depth_] + // Above
                    pressure_[index - width_] + // Behind
                    pressure_[index + width_]    // In front
                    ) / 6.0f;
            }
        }
    }
    pressure_ = newPressure_; // Update pressure
}


void CFDWaveBox::calculateVelocity(float deltaTime) {
    // Simple velocity calculation (example: proportional to pressure gradient)
    for (int x = 1; x < width_ - 1; ++x) {
        for (int y = 1; y < height_ - 1; ++y) {
            for (int z = 1; z < depth_ - 1; ++z) {
                int index = x + y * width_ * depth_ + z * width_;

                // Calculate pressure gradients
                float dpdx = (pressure_[index + 1] - pressure_[index - 1]) / (2 * gridSize_);
                float dpdy = (pressure_[index + width_ * depth_] - pressure_[index - width_ * depth_]) / (2 * gridSize_);
                float dpdz = (pressure_[index + width_] - pressure_[index - width_]) / (2 * gridSize_);

                // Update velocities (very simplified Navier-Stokes)
                velocityX_[index] -= dpdx * deltaTime / density_;
                velocityY_[index] -= dpdy * deltaTime / density_;
                velocityZ_[index] -= dpdz * deltaTime / density_;

                // Apply viscosity (simple damping)
                velocityX_[index] *= (1.0f - viscosity_ * deltaTime);
                velocityY_[index] *= (1.0f - viscosity_ * deltaTime);
                velocityZ_[index] *= (1.0f - viscosity_ * deltaTime);
            }
        }
    }
}

*/