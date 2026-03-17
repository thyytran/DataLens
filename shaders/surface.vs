#version 330 core

/*
 * surface.vs - Molecular Surface Vertex Shader
 * 
 * Input: Surface mesh vertices with position, normal, and color
 * Output: Transformed position and data for surface shading
 */

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec4 aColor;

out vec3 FragPos;
out vec3 Normal;
out vec4 Color;
out float Depth;

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;

void main() {
    // Transform position to world space
    vec4 worldPos = u_model * vec4(aPos, 1.0);
    FragPos = worldPos.xyz;
    
    // Transform normal to world space
    mat3 normalMatrix = transpose(inverse(mat3(u_model)));
    Normal = normalize(normalMatrix * aNormal);
    
    // Pass through color
    Color = aColor;
    
    // Clip-space position
    vec4 clipPos = u_projection * u_view * worldPos;
    gl_Position = clipPos;
    
    // Normalized depth for depth-based effects
    Depth = clipPos.z / clipPos.w * 0.5 + 0.5;
}
