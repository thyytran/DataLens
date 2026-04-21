#version 330 core

/*
 * ribbon.fs - Ribbon/Cartoon Representation Fragment Shader
 * 
 * Features:
 * - Two-sided lighting (for flat ribbons)
 * - Blinn-Phong shading with configurable components
 * - Rim/fresnel lighting for edge highlighting
 * - Per-vertex color from secondary structure
 */

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in vec4 Color;

out vec4 FragColor;

// Light and material uniforms
uniform vec3 u_lightDir;        // Normalized light direction (world space)
uniform vec3 u_cameraPos;       // Camera position (world space)
uniform vec3 u_lightColor;      // Light color (default: white)
uniform float u_ambientStrength;    // Ambient coefficient (default: 0.15)
uniform float u_diffuseStrength;    // Diffuse coefficient (default: 0.7)
uniform float u_specularStrength;   // Specular coefficient (default: 0.4)
uniform float u_shininess;          // Specular exponent (default: 32.0)
uniform float u_rimStrength;        // Rim light strength (default: 0.2)
uniform float u_rimPower;           // Rim light falloff (default: 3.0)

void main() {
    // Normalize input vectors
    vec3 normal = normalize(Normal);
    vec3 lightDir = normalize(u_lightDir);
    vec3 viewDir = normalize(u_cameraPos - FragPos);
    
    // Two-sided lighting: flip normal if facing away from camera
    if (dot(normal, viewDir) < 0.0) {
        normal = -normal;
    }
    
    // Ambient component
    float ambient = u_ambientStrength;
    
    // Diffuse component (Lambert)
    float diff = max(dot(normal, lightDir), 0.0);
    float diffuse = u_diffuseStrength * diff;
    
    // Specular component (Blinn-Phong)
    vec3 halfDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfDir), 0.0), u_shininess);
    float specular = u_specularStrength * spec;
    
    // Rim lighting (fresnel-like edge glow)
    float rim = 1.0 - max(dot(normal, viewDir), 0.0);
    rim = pow(rim, u_rimPower) * u_rimStrength;
    
    // Combine lighting components
    vec3 lighting = vec3(ambient + diffuse + rim) * Color.rgb + 
                    vec3(specular) * u_lightColor;
    
    // Add slight variation based on UV for visual interest
    float uvVariation = 1.0 + 0.05 * sin(TexCoord.y * 10.0);
    
    FragColor = vec4(lighting * uvVariation, Color.a);
    
    // Gamma correction (if not using sRGB framebuffer)
    // FragColor.rgb = pow(FragColor.rgb, vec3(1.0/2.2));
}
