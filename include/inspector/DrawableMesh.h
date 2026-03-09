#pragma once

#include <../lib/glad/glad.h>
#include <GLFW/glfw3.h>
#include "../lib/obj/OBJ_Loader.h"

class DrawableMesh
{
    unsigned int VAO;
    unsigned int VBO;
    unsigned int EBO;
    unsigned int texture;

public:
    unsigned int vert_count;
    unsigned int ind_count;

    DrawableMesh(GLuint drawMode,
        float* vertices, unsigned int vert_count,
        unsigned int* indices, unsigned int ind_count, const char* texture_path)
        : DrawableMesh(drawMode, vertices, vert_count, indices, ind_count, false, texture_path) {
    };

    DrawableMesh(GLuint draw_mode, float* vertices, unsigned int vertex_count,
        unsigned int* indices, unsigned int indices_count,
        bool color = false, const char* texture_path = nullptr);

    DrawableMesh(GLuint drawMode, objl::Mesh mesh, const char* texturesFolder = nullptr);

    void Draw() const;
    void LoadTexture(const char* texture_path);
};

