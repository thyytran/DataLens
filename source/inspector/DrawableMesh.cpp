#define STB_IMAGE_IMPLEMENTATION

#include <iostream>
#include <../lib/glad/glad.h>
#include "../include/graphics/stb_image.h"
#include "../include/inspector/DrawableMesh.h"
#include "../lib/obj/OBJ_Loader.h"

/**
 * @brief Construct a new Mesh object with vertex and index data.
*/
DrawableMesh::DrawableMesh(GLuint draw_mode,
    float* vertices, unsigned int vertex_count,
    unsigned int* indices, unsigned int indices_count,
    bool color, const char* texture_path) {
    // Stores number of vertices and indices in the object to member variables named vert_count and ind_count
    this->vert_count = vertex_count;
    this->ind_count = indices_count;

    // Generate Vertex Array Object
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    // Generate Vertex Buffer Object
    glGenBuffers(1, &VBO);

    // Binds the above VBO to GL_ARRAY_BUFFER target. Making the VBO the active buffer for that target
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    // Allocates memory for the VBO and copies vertex data into it
    glBufferData(GL_ARRAY_BUFFER, vertex_count * sizeof(float), vertices, draw_mode);

    // Generates a new buffer object and stores its ID in EBO
    glGenBuffers(1, &EBO);

    // Binds the newly created buffer to GL_ELEMENT_ARRAY_BUFFER target, indicating for index data
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);

    // Allocates storage and uploads the index data to the GPU
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices_count * sizeof(float), indices, draw_mode);

    // Vertex Attribute Layout (VAL) - this defines constants for vertex attribute locations in the shader
    const unsigned int v_attribute = 0; // For vertex positions
    const unsigned int c_attribute = 1; // For vertex colors
    const unsigned int t_attribute = 2; // For texture coordinates
    const unsigned int n_attribute = 3; // For normal vectors

    if (texture_path != nullptr) {
        LoadTexture(texture_path);
    }

    // Stride calculation - initially set to 3 (for xyz position)
    int vertex_color_stride = 3;

    if (color) {
        // If color data is present, the stride is increased to 7
        vertex_color_stride = 7;

        // If texture data is present, stride increases by 2. Sets up texture coordinate attribute pointer
        if (texture_path != nullptr) {
            vertex_color_stride += 2;
            glVertexAttribPointer(t_attribute, 2,
                GL_FLOAT, GL_FALSE,
                vertex_color_stride * sizeof(float), (void*)(7 * sizeof(float)));
            glEnableVertexAttribArray(t_attribute);
        }

        // Set up color attribute pointers
        glVertexAttribPointer(c_attribute, 4,
            GL_FLOAT, GL_FALSE,
            vertex_color_stride * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(c_attribute);
    }
    else {
        // If no color but texture is present, sets up texture coordinate attribute pointer differently
        if (texture_path != nullptr) {
            // move 5 each times stride, start at 3
            glVertexAttribPointer(t_attribute, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(t_attribute);
        }
    }
    // Setup position attribute pointer always present
    glVertexAttribPointer(v_attribute, 3,
        GL_FLOAT, GL_FALSE,
        vertex_color_stride * sizeof(float), (void*)0);
    glEnableVertexAttribArray(v_attribute);
}

/*
 * Constructor for a DrawableMesh class
 * Renders 3D meshes loaded using the objl lib (a Wavefront OBJ file loader)
 */
DrawableMesh::DrawableMesh(GLuint drawMode, objl::Mesh mesh, const char* texturesFolder) {

    // Store the number of vertices and indices from the input objl::Mesh obj
    // into member variable of DrawableMesh class
    this->vert_count = mesh.Vertices.size();
    this->ind_count = mesh.Indices.size();

    // Generates and binds a Vertex Array Object
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    // Generates and binds a Vertex Buffer Object
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO); // bind vertex buffer

    glBufferData(GL_ARRAY_BUFFER, vert_count * sizeof(objl::Vertex), &mesh.Vertices[0], drawMode);
    // send vertex data to buffer

    // Generates an Element Buffer Object (EBO), bind it, and upload the index data to the EBO
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO); // bind index buffer
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, ind_count * sizeof(unsigned int), &mesh.Indices[0], drawMode);

    // Defines constants for the Vertex Attribute Locations
    const unsigned int v_attribute = 0; // Position
    const unsigned int t_attribute = 1; // Texture coordinates
    const unsigned int n_attribute = 2; // Normals

    // Texture loading logic
    // Load a texture based on the mesh.MeshMaterial.map_Kd value (diffuse texture file name)
    auto file = mesh.MeshMaterial.map_Kd;

    //std::cout << "Mesh: " << mesh.MeshName << " - material texture: " << file << std::endl;

    // Checks if a texture folder path was provided
    if (texturesFolder != nullptr) {

        // Checks if a texture file name is available (presumably from the .mtl file)
        if (!file.empty()) {
            auto const true_path = std::string(texturesFolder) + file;
            std::vector<char> charArray(true_path.size() + 1);
            // std::strcpy(charArray.data(), true_path.c_str()); // Use charArray.data() to get a char* // uncomment when fixing strcpy
            std::cout << "Path: " << charArray.data() << std::endl;
            LoadTexture(charArray.data());
        }
        else {
            std::cout << "No texture specified" << std::endl;
            LoadTexture("resources/black.png");
        }
    }
    else {
        std::cout << "No texture folder specified" << std::endl;
        LoadTexture("resources/white.png");
    }

    // v_attribute: this attribute holds the 3D position (x,y,z coordinates) of each vertex
    // Defining shape of the model
    glVertexAttribPointer(v_attribute, 3,
        GL_FLOAT, GL_FALSE, sizeof(objl::Vertex), (void*)0);
    glEnableVertexAttribArray(v_attribute);

    // Tells OpenGL how to interpret the vertex data within the VBO
    glVertexAttribPointer(t_attribute, 2,
        GL_FLOAT, GL_FALSE,
        sizeof(objl::Vertex), (void*)offsetof(objl::Vertex, TextureCoordinate));

    // Enables the specified vertex attribute so that it's used during rendering
    // t_attribute: this attribute holds the 2D texture coordinates (u, v values) for each vertex
    // Texture coordinates map the vertices to points on a 2D texture image
    glEnableVertexAttribArray(t_attribute);

    // n_attribute: this attribute stores normal vector for each vertex, crucial for lighting calculations
    // Defining the direction a surface is facing, which determines how light interacts with it.
    glVertexAttribPointer(n_attribute, 3,
        GL_FLOAT, GL_FALSE,
        sizeof(objl::Vertex), (void*)offsetof(objl::Vertex, Normal));
    glEnableVertexAttribArray(n_attribute);
}

/*
 *  Defines the Draw method for class DrawableMesh
 *  Draws a mesh using OpenGL
 */
void DrawableMesh::Draw() const {
    // Binds a texture to GL_TEXTURE_2D -> subsequent drawing operation will use this texture
    glBindTexture(GL_TEXTURE_2D, texture);

    // Binds the Vertex Array Object associated with the mesh
    glBindVertexArray(VAO);

    // Core drawing command, draw the mesh using currently bound VAO and specified indices
    // nullptr is the offset into the index buffer, nullptr or 0 means to start from the beginning of the buffer
    glDrawElements(GL_TRIANGLES, ind_count, GL_UNSIGNED_INT, nullptr);
}

/*
 * Defines the LoadTexture method for class DrawableMesh
 * Loads and sets up a texture from a file and sets it up for use in OpenGL using sbti_image
 */
void DrawableMesh::LoadTexture(const char* texture_path) {
    // Image loading, images loaded from files are typically stored with the y-axis flipped
    // Tells stb_image lib to flip the image vertically during loading
    stbi_set_flip_vertically_on_load(true);

    // Generate a new texture object, creates 1 or more textures names (IDs)
    // texture should be a GLuint member of the DrawableMesh class
    glGenTextures(1, &texture);

    // Binds the newly generated texture to GL_TEXTURE_2D target
    glBindTexture(GL_TEXTURE_2D, texture);

    // Set the texture wrapping mode for the s and t axis (x and y axes in image space)
    // GL_REPEAT means the texture will repeat if the texture coordinates go outside the range [0, 1]
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    // Texture filtering parameters
    // GL_TEXTURE_MIN_FILTER controls how the texture is sampled when it's minified
    // GL_TEXTURE_MAG_FILTER controls how the texture is sampled when it's magnified
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    int width, height, numberOfChannels;

    // 0 as last arg tells sbti_load to load the image with any number of channels
    unsigned char* texture_data = stbi_load(texture_path, &width, &height, &numberOfChannels, 0);
    if (texture_data) {
        std::cout << "Texture loaded!" << std::endl;

        // Sets the pixel storage alignment
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

        // Uploads texture data to the GPU
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height,
            0, GL_RGB, GL_UNSIGNED_BYTE, texture_data);

        // Generate mipmaps for the textures
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    else {
        // Error handling
        std::cout << "Texture loading failed. " << std::endl;
    }
    stbi_image_free(texture_data);
}
