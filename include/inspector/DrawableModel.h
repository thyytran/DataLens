#pragma once

#ifndef OPENGL_MODEL_VIEWER_DRAWABLE_MODEL_H
#define OPENGL_MODEL_VIEWER_DRAWABLE_MODEL_H

#include <../lib/glad/glad.h>
#include <GLFW/glfw3.h>
#include "../include/inspector/DrawableMesh.h"
#include "../lib/glm/vec3.hpp"

/**
 * Served as a class to represent a 3D model in
 */
class DrawableModel {
public:
    // Constructors and methods
    DrawableModel(GLuint drawMode, const char* objPath, const char* texturesFolder = nullptr);
    void Draw();

    // Class members
    Vec3 avg_pos;// A 3D vector representing the average position of the model
    // glm::vec3 avg_pos;
    std::vector<DrawableMesh> meshes; // A vector of DrawableMesh objects, representing individual parts or components
    unsigned int mesh_count; // Number of meshes in the model
    unsigned int vertex_count; // Number of vertices in the model
    unsigned int material_count; // Number of material used in the models
};

#endif //OPENGL_MODEL_VIEWER_DRAWABLE_MODEL_H
