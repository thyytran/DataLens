#pragma once

#include <vector>
#include <memory>
#include <set>
#include "RepresentationType.h"
#include "RibbonTemplate.h"
#include "SurfaceTemplate.h"
#include "bio/SecondaryStructure.h"

// Forward declarations (match your existing project structure)
class Shader;
class Window;
class Camera;
class MoleculeData;

/*
 * RepresentationRenderer.h
 * Unified renderer for all molecular representations.
 * 
 * Manages:
 * - Ball-and-Stick (delegates to existing Model class)
 * - Ribbon/Cartoon (RibbonTemplate)
 * - Surface (SurfaceTemplate)
 * 
 * Usage in render loop:
 *   renderer.render(shaders, window, camera, modelMatrix);
 */

class RepresentationRenderer {
public:
    // Rendering parameters
    float ballRadius = 0.4f;          // Radius for ball representation
    float stickRadius = 0.15f;        // Radius for stick/bond representation
    float surfaceOpacity = 0.8f;      // Surface transparency
    
private:
    // Representation generators
    RibbonTemplate ribbonTemplate;
    SurfaceTemplate surfaceTemplate;
    
    // Active representations (bit flags)
    RepresentationFlags activeRepresentations = RepresentationFlags::BALL_AND_STICK;
    
    // Cached molecule data pointer
    MoleculeData* cachedMoleculeData = nullptr;
    
    // OpenGL resources for ribbon
    unsigned int ribbonVAO = 0;
    unsigned int ribbonVBO = 0;
    unsigned int ribbonEBO = 0;
    bool ribbonBuffersValid = false;
    
    // OpenGL resources for surface
    unsigned int surfaceVAO = 0;
    unsigned int surfaceVBO = 0;
    unsigned int surfaceEBO = 0;
    bool surfaceBuffersValid = false;
    
    // Surface type
    SurfaceTemplate::SurfaceType surfaceType = SurfaceTemplate::SurfaceType::SOLVENT_ACCESSIBLE;
    
public:
    RepresentationRenderer();
    ~RepresentationRenderer();
    
    // Non-copyable
    RepresentationRenderer(const RepresentationRenderer&) = delete;
    RepresentationRenderer& operator=(const RepresentationRenderer&) = delete;
    
    /*
     * Initialize OpenGL resources.
     * Call after OpenGL context is created.
     */
    void initialize();
    
    /*
     * Clean up OpenGL resources.
     * Call before destroying OpenGL context.
     */
    void cleanup();
    
    /*
     * Load molecule data and generate all representations.
     * @param moleculeData Pointer to loaded molecule (PDBFile, etc.)
     */
    void loadMoleculeData(MoleculeData* moleculeData);
    
    /*
     * Enable or disable a representation type.
     */
    void setRepresentationEnabled(RepresentationType type, bool enabled);
    
    /*
     * Check if a representation is enabled.
     */
    bool isRepresentationEnabled(RepresentationType type) const;
    
    /*
     * Toggle a representation on/off.
     */
    void toggleRepresentation(RepresentationType type);
    
    /*
     * Set active representations using flags.
     */
    void setActiveRepresentations(RepresentationFlags flags);
    
    /*
     * Get current active representation flags.
     */
    RepresentationFlags getActiveRepresentations() const;
    
    /*
     * Set surface type (VDW, SAS, SES).
     */
    void setSurfaceType(SurfaceTemplate::SurfaceType type);
    
    /*
     * Regenerate ribbon geometry (e.g., after secondary structure update).
     */
    void regenerateRibbon();
    
    /*
     * Regenerate surface geometry (e.g., after atom selection change).
     */
    void regenerateSurface();
    
    /*
     * Render all active representations.
     * @param ribbonShader Shader for ribbon rendering
     * @param surfaceShader Shader for surface rendering
     * @param window Window for viewport info
     * @param camera Camera for view/projection matrices
     * @param modelMatrix Model transformation matrix (from Model::getModelMatrix)
     */
    void render(
        Shader* ribbonShader,
        Shader* surfaceShader,
        const Window* window,
        const Camera* camera,
        const float* modelMatrix  // 4x4 matrix as float[16]
    );
    
    /*
     * Get the number of vertices in the ribbon mesh.
     */
    size_t getRibbonVertexCount() const;
    
    /*
     * Get the number of vertices in the surface mesh.
     */
    size_t getSurfaceVertexCount() const;
    
private:
    /*
     * Extract ribbon control points from chains.
     */
    void extractRibbonControlPoints(std::vector<RibbonControlPoint>& controlPoints);
    
    /*
     * Extract surface atoms from all atoms.
     */
    void extractSurfaceAtoms(std::vector<SurfaceAtom>& surfaceAtoms);
    
    /*
     * Upload ribbon geometry to GPU.
     */
    void uploadRibbonBuffers();
    
    /*
     * Upload surface geometry to GPU.
     */
    void uploadSurfaceBuffers();
    
    /*
     * Render ribbon representation.
     */
    void renderRibbon(Shader* shader, const Window* window, const Camera* camera, 
                      const float* modelMatrix);
    
    /*
     * Render surface representation.
     */
    void renderSurface(Shader* shader, const Window* window, const Camera* camera,
                       const float* modelMatrix);
    
    /*
     * Get secondary structure color.
     */
    void getSSColor(SecondaryStructure ss, float& r, float& g, float& b);
    
    /*
     * Get secondary structure type for a residue using helix/sheet data.
     */
    SecondaryStructure getSecondaryStructureForResidue(char chainId, int residueNum) const;
};
