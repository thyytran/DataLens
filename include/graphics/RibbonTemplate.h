#pragma once

#include <vector>
#include "../bio/SecondaryStructure.h"

/*
 * RibbonTemplate.h
 * Generates ribbon/cartoon geometry for protein secondary structure visualization.
 * 
 * The ribbon representation shows:
 * - Helices as wide, flat ribbons (coiled)
 * - Sheets as flat arrows pointing in strand direction
 * - Coils/loops as thin tubes
 * 
 * Uses Catmull-Rom splines for smooth interpolation between C-alpha positions.
 */

// Control point for ribbon spline (extracted from residue C-alphas)
struct RibbonControlPoint {
    float x, y, z;           // Position (C-alpha)
    float nx, ny, nz;        // Guide normal (C-alpha to C-beta direction)
    SecondaryStructure ssType;
    int residueIndex;
    char chainId;
    float colorR, colorG, colorB, colorA;
};

// Vertex for ribbon mesh
struct RibbonVertex {
    float x, y, z;           // Position
    float nx, ny, nz;        // Normal
    float u, v;              // Texture coordinates
    float r, g, b, a;        // Color
};

class RibbonTemplate {
public:
    // Ribbon geometry parameters
    static constexpr float HELIX_WIDTH = 2.0f;       // Width of helix ribbon
    static constexpr float HELIX_THICKNESS = 0.4f;   // Thickness of helix ribbon
    static constexpr float SHEET_WIDTH = 2.2f;       // Width of sheet ribbon
    static constexpr float SHEET_THICKNESS = 0.5f;   // Thickness of sheet ribbon
    static constexpr float SHEET_ARROW_WIDTH = 3.0f; // Width at arrow head
    static constexpr float COIL_RADIUS = 0.3f;       // Radius of coil tube
    static constexpr int SPLINE_SEGMENTS = 4;        // Segments per residue
    static constexpr int TUBE_SIDES = 8;             // Sides for coil tube
    
private:
    std::vector<RibbonVertex> vertices;
    std::vector<unsigned int> indices;
    
public:
    RibbonTemplate() = default;
    
    /*
     * Generate ribbon geometry from control points.
     * Control points should be extracted from C-alpha positions along the backbone.
     */
    void generateRibbon(const std::vector<RibbonControlPoint>& controlPoints);
    
    /*
     * Clear all generated geometry.
     */
    void clear();
    
    /*
     * Get generated vertices.
     */
    const std::vector<RibbonVertex>& getVertices() const { return vertices; }
    
    /*
     * Get generated indices.
     */
    const std::vector<unsigned int>& getIndices() const { return indices; }
    
    /*
     * Get vertex count.
     */
    size_t vertexCount() const { return vertices.size(); }
    
    /*
     * Get index count.
     */
    size_t indexCount() const { return indices.size(); }
    
    /*
     * Get vertex data pointer for OpenGL upload.
     */
    const float* getVertexData() const {
        return vertices.empty() ? nullptr : &vertices[0].x;
    }
    
    /*
     * Get index data pointer for OpenGL upload.
     */
    const unsigned int* getIndexData() const {
        return indices.empty() ? nullptr : indices.data();
    }
    
    /*
     * Vertex stride in bytes.
     */
    static constexpr size_t VERTEX_STRIDE = sizeof(RibbonVertex);
    
private:
    /*
     * Catmull-Rom spline interpolation.
     * Returns position at parameter t [0,1] between points p1 and p2.
     */
    void catmullRom(
        float p0x, float p0y, float p0z,
        float p1x, float p1y, float p1z,
        float p2x, float p2y, float p2z,
        float p3x, float p3y, float p3z,
        float t,
        float& outX, float& outY, float& outZ
    );
    
    /*
     * Catmull-Rom tangent (derivative).
     */
    void catmullRomTangent(
        float p0x, float p0y, float p0z,
        float p1x, float p1y, float p1z,
        float p2x, float p2y, float p2z,
        float p3x, float p3y, float p3z,
        float t,
        float& outX, float& outY, float& outZ
    );
    
    /*
     * Generate cross-section vertices for helix (flat ribbon).
     */
    void generateHelixCrossSection(
        float cx, float cy, float cz,
        float tx, float ty, float tz,      // Tangent
        float nx, float ny, float nz,      // Normal (up direction)
        float width, float thickness,
        float r, float g, float b, float a,
        float v,
        std::vector<RibbonVertex>& crossSection
    );
    
    /*
     * Generate cross-section vertices for sheet (flat ribbon with arrow).
     */
    void generateSheetCrossSection(
        float cx, float cy, float cz,
        float tx, float ty, float tz,
        float nx, float ny, float nz,
        float width, float thickness,
        float r, float g, float b, float a,
        float v,
        bool isArrowHead,
        std::vector<RibbonVertex>& crossSection
    );
    
    /*
     * Generate cross-section vertices for coil (tube).
     */
    void generateCoilCrossSection(
        float cx, float cy, float cz,
        float tx, float ty, float tz,
        float nx, float ny, float nz,
        float radius,
        float r, float g, float b, float a,
        float v,
        std::vector<RibbonVertex>& crossSection
    );
    
    /*
     * Connect two cross-sections with triangles.
     */
    void connectCrossSections(
        const std::vector<RibbonVertex>& cs1,
        const std::vector<RibbonVertex>& cs2,
        unsigned int baseIndex1,
        unsigned int baseIndex2
    );
    
    /*
     * Normalize a vector in place.
     */
    void normalize(float& x, float& y, float& z);
    
    /*
     * Cross product.
     */
    void cross(
        float ax, float ay, float az,
        float bx, float by, float bz,
        float& rx, float& ry, float& rz
    );
};
