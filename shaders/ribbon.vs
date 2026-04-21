#version 330 core

/*
 * ribbon.vs - Ribbon/Cartoon Representation Vertex Shader
 * 
 * Input: Ribbon mesh vertices with position, normal, UV, and color
 * Output: Transformed position and data for Blinn-Phong shading
 */

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;
layout (location = 3) in vec4 aColor;

out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;
out vec4 Color;

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;

void main() {
    // Transform position to world space
    vec4 worldPos = u_model * vec4(aPos, 1.0);
    FragPos = worldPos.xyz;
    
    // Transform normal to world space (use normal matrix for non-uniform scaling)
    mat3 normalMatrix = transpose(inverse(mat3(u_model)));
    Normal = normalize(normalMatrix * aNormal);
    
    // Pass through texture coordinates and color
    TexCoord = aTexCoord;
    Color = aColor;
    
    // Final clip-space position
    gl_Position = u_projection * u_view * worldPos;
}
