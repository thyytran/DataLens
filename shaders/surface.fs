#version 330 core

/*
 * surface.fs - Molecular Surface Fragment Shader
 * 
 * Features:
 * - Semi-transparent rendering
 * - Fresnel effect (edges more visible)
 * - Wrap lighting for soft appearance
 * - Depth-based opacity for interior clarity
 */

in vec3 FragPos;
in vec3 Normal;
in vec4 Color;
in float Depth;

out vec4 FragColor;

// Uniforms
uniform vec3 u_lightDir;
uniform vec3 u_cameraPos;
uniform vec3 u_lightColor;
uniform float u_opacity;            // Base opacity (default: 0.8)
uniform float u_ambientStrength;    // Ambient (default: 0.2)
uniform float u_diffuseStrength;    // Diffuse (default: 0.6)
uniform float u_specularStrength;   // Specular (default: 0.3)
uniform float u_shininess;          // Shininess (default: 16.0)
uniform float u_fresnelPower;       // Fresnel falloff (default: 2.0)
uniform float u_fresnelStrength;    // Fresnel intensity (default: 0.4)
uniform float u_wrapAmount;         // Wrap lighting amount (default: 0.3)

void main() {
    // Normalize vectors
    vec3 normal = normalize(Normal);
    vec3 lightDir = normalize(u_lightDir);
    vec3 viewDir = normalize(u_cameraPos - FragPos);
    
    // Two-sided normals
    if (dot(normal, viewDir) < 0.0) {
        normal = -normal;
    }
    
    // Ambient
    float ambient = u_ambientStrength;
    
    // Wrap diffuse lighting (softer shadows)
    float NdotL = dot(normal, lightDir);
    float diffuse = u_diffuseStrength * max((NdotL + u_wrapAmount) / (1.0 + u_wrapAmount), 0.0);
    
    // Specular (Blinn-Phong)
    vec3 halfDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(normal, halfDir), 0.0), u_shininess);
    float specular = u_specularStrength * spec;
    
    // Fresnel effect (Schlick approximation)
    // Makes edges more opaque/reflective
    float fresnel = 1.0 - max(dot(normal, viewDir), 0.0);
    fresnel = pow(fresnel, u_fresnelPower) * u_fresnelStrength;
    
    // Combine lighting
    vec3 baseColor = Color.rgb;
    vec3 litColor = baseColor * (ambient + diffuse) + u_lightColor * specular;
    
    // Add fresnel rim effect
    litColor += fresnel * u_lightColor * 0.5;
    
    // Compute alpha with fresnel influence
    // Edges are more opaque due to Fresnel
    float alpha = u_opacity + fresnel * (1.0 - u_opacity) * 0.5;
    alpha = min(alpha, 1.0);
    
    // Slight depth-based transparency (closer to camera = more opaque)
    // This helps see interior structure
    // alpha *= mix(0.8, 1.0, 1.0 - Depth);
    
    FragColor = vec4(litColor, alpha);
    
    // Premultiplied alpha for better blending
    FragColor.rgb *= FragColor.a;
}
