#pragma once

#include <vector>

/*
 * SurfaceTemplate.h
 * Generates molecular surface geometry using the Marching Cubes algorithm.
 * 
 * Supports three surface types:
 * - Van der Waals (VDW): Union of atomic spheres
 * - Solvent Accessible Surface (SAS): VDW expanded by probe radius
 * - Solvent Excluded Surface (SES): The molecular surface (Connolly surface)
 */

// Input atom for surface generation
struct SurfaceAtom {
    float x, y, z;           // Position
    float radius;            // Van der Waals radius
    float r, g, b, a;        // Color
    int residueIndex;
    char chainId;
};

// Output vertex for surface mesh
struct SurfaceVertex {
    float x, y, z;           // Position
    float nx, ny, nz;        // Normal
    float r, g, b, a;        // Color
};

class SurfaceTemplate {
public:
    // Surface type enumeration
    enum class SurfaceType {
        VAN_DER_WAALS,       // Union of atomic spheres
        SOLVENT_ACCESSIBLE,  // VDW + probe radius
        SOLVENT_EXCLUDED     // Connolly molecular surface
    };
    
    // Default parameters
    static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;    // Water probe radius (Angstroms)
    static constexpr float DEFAULT_GRID_SPACING = 0.5f;    // Grid resolution
    static constexpr float DEFAULT_ISOVALUE = 0.0f;        // Surface isovalue
    
private:
    std::vector<SurfaceVertex> vertices;
    std::vector<unsigned int> indices;
    
    // Grid parameters
    float gridSpacing = DEFAULT_GRID_SPACING;
    float probeRadius = DEFAULT_PROBE_RADIUS;
    float isovalue = DEFAULT_ISOVALUE;
    
    // Grid dimensions
    int gridX = 0, gridY = 0, gridZ = 0;
    float minX, minY, minZ;
    float maxX, maxY, maxZ;
    
    // Scalar field and color field
    std::vector<float> scalarField;
    std::vector<float> colorField;  // RGBA interleaved
    
public:
    SurfaceTemplate() = default;
    
    /*
     * Generate molecular surface from atoms.
     */
    void generateSurface(const std::vector<SurfaceAtom>& atoms,
                         SurfaceType type = SurfaceType::SOLVENT_ACCESSIBLE);
    
    /*
     * Set grid spacing (smaller = higher resolution, more memory).
     */
    void setGridSpacing(float spacing) { gridSpacing = spacing; }
    
    /*
     * Set probe radius for SAS/SES surfaces.
     */
    void setProbeRadius(float radius) { probeRadius = radius; }
    
    /*
     * Clear all generated geometry.
     */
    void clear();
    
    /*
     * Get generated vertices.
     */
    const std::vector<SurfaceVertex>& getVertices() const { return vertices; }
    
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
    static constexpr size_t VERTEX_STRIDE = sizeof(SurfaceVertex);
    
private:
    /*
     * Compute bounding box from atoms.
     */
    void computeBoundingBox(const std::vector<SurfaceAtom>& atoms, float padding);
    
    /*
     * Initialize the 3D scalar field grid.
     */
    void initializeGrid();
    
    /*
     * Compute scalar field values for each grid point.
     */
    void computeScalarField(const std::vector<SurfaceAtom>& atoms, SurfaceType type);
    
    /*
     * Run Marching Cubes to extract the isosurface.
     */
    void marchingCubes();
    
    /*
     * Process a single cube in the Marching Cubes algorithm.
     */
    void processCell(int ix, int iy, int iz);
    
    /*
     * Get scalar field value at grid index.
     */
    float getFieldValue(int ix, int iy, int iz) const;
    
    /*
     * Get color at grid index.
     */
    void getFieldColor(int ix, int iy, int iz, float& r, float& g, float& b, float& a) const;
    
    /*
     * Interpolate vertex position along an edge.
     */
    void interpolateVertex(
        float x1, float y1, float z1, float v1,
        float x2, float y2, float z2, float v2,
        float& outX, float& outY, float& outZ
    );
    
    /*
     * Interpolate color along an edge.
     */
    void interpolateColor(
        float r1, float g1, float b1, float a1, float v1,
        float r2, float g2, float b2, float a2, float v2,
        float& outR, float& outG, float& outB, float& outA
    );
    
    /*
     * Compute gradient (normal) at a point using central differences.
     */
    void computeGradient(float x, float y, float z, float& nx, float& ny, float& nz);
    
    /*
     * Normalize a vector.
     */
    void normalize(float& x, float& y, float& z);
    
    // Marching Cubes lookup tables
    static const int edgeTable[256];
    static const int triTable[256][16];
};
