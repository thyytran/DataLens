#include "graphics/RepresentationRenderer.h"
#include <GL/glew.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>

// Include actual project headers
#include "bio/MoleculeData.h"
#include "bio/Helix.h"
#include "bio/Sheet.h"

/*
 * RepresentationRenderer.cpp
 * Implementation of the unified molecular representation renderer.
 * Works with the actual MoleculeData structure from DataLens.
 */

// VDW radii for common elements (Angstroms)
static float getVDWRadius(const std::string& element) {
    if (element == "C" || element == "c") return 1.70f;
    if (element == "N" || element == "n") return 1.55f;
    if (element == "O" || element == "o") return 1.52f;
    if (element == "S" || element == "s") return 1.80f;
    if (element == "H" || element == "h") return 1.20f;
    if (element == "P" || element == "p") return 1.80f;
    if (element == "FE" || element == "Fe" || element == "fe") return 1.95f;
    if (element == "ZN" || element == "Zn" || element == "zn") return 1.39f;
    if (element == "CA" || element == "Ca" || element == "ca") return 1.97f;
    if (element == "MG" || element == "Mg" || element == "mg") return 1.73f;
    return 1.50f;  // Default
}

// CPK coloring scheme
static void getCPKColor(const std::string& element, float& r, float& g, float& b) {
    if (element == "C" || element == "c") { r = 0.5f; g = 0.5f; b = 0.5f; return; }
    if (element == "N" || element == "n") { r = 0.2f; g = 0.2f; b = 0.8f; return; }
    if (element == "O" || element == "o") { r = 0.8f; g = 0.2f; b = 0.2f; return; }
    if (element == "S" || element == "s") { r = 0.8f; g = 0.8f; b = 0.2f; return; }
    if (element == "H" || element == "h") { r = 0.9f; g = 0.9f; b = 0.9f; return; }
    if (element == "P" || element == "p") { r = 0.8f; g = 0.5f; b = 0.0f; return; }
    r = 0.7f; g = 0.7f; b = 0.7f;  // Default gray
}

RepresentationRenderer::RepresentationRenderer() 
    : cachedMoleculeData(nullptr) {}

RepresentationRenderer::~RepresentationRenderer() {
    cleanup();
}

void RepresentationRenderer::initialize() {
    // Generate ribbon buffers
    glGenVertexArrays(1, &ribbonVAO);
    glGenBuffers(1, &ribbonVBO);
    glGenBuffers(1, &ribbonEBO);
    
    // Generate surface buffers
    glGenVertexArrays(1, &surfaceVAO);
    glGenBuffers(1, &surfaceVBO);
    glGenBuffers(1, &surfaceEBO);
}

void RepresentationRenderer::cleanup() {
    if (ribbonVAO) {
        glDeleteVertexArrays(1, &ribbonVAO);
        ribbonVAO = 0;
    }
    if (ribbonVBO) {
        glDeleteBuffers(1, &ribbonVBO);
        ribbonVBO = 0;
    }
    if (ribbonEBO) {
        glDeleteBuffers(1, &ribbonEBO);
        ribbonEBO = 0;
    }
    
    if (surfaceVAO) {
        glDeleteVertexArrays(1, &surfaceVAO);
        surfaceVAO = 0;
    }
    if (surfaceVBO) {
        glDeleteBuffers(1, &surfaceVBO);
        surfaceVBO = 0;
    }
    if (surfaceEBO) {
        glDeleteBuffers(1, &surfaceEBO);
        surfaceEBO = 0;
    }
    
    ribbonBuffersValid = false;
    surfaceBuffersValid = false;
}

void RepresentationRenderer::getSSColor(SecondaryStructure ss, float& r, float& g, float& b) {
    switch (ss) {
        case SecondaryStructure::HELIX_ALPHA:
        case SecondaryStructure::HELIX_310:
        case SecondaryStructure::HELIX_PI:
            r = 1.0f; g = 0.2f; b = 0.2f;  // Red for helices
            break;
        case SecondaryStructure::SHEET:
        case SecondaryStructure::BRIDGE:
            r = 1.0f; g = 1.0f; b = 0.2f;  // Yellow for sheets
            break;
        case SecondaryStructure::TURN:
            r = 0.2f; g = 0.8f; b = 0.8f;  // Cyan for turns
            break;
        case SecondaryStructure::COIL:
        default:
            r = 0.7f; g = 0.7f; b = 0.7f;  // Gray for coil
            break;
    }
}

SecondaryStructure RepresentationRenderer::getSecondaryStructureForResidue(
    char chainId, int residueNum
) const {
    if (!cachedMoleculeData) return SecondaryStructure::COIL;
    
    // Check helices
    for (const Helix& helix : cachedMoleculeData->helices) {
        if (helix.chain == chainId && 
            residueNum >= helix.residueStart && 
            residueNum <= helix.residueEnd) {
            // Convert helix type to SecondaryStructure
            switch (helix.type) {
                case 1: return SecondaryStructure::HELIX_ALPHA;  // Right-handed alpha
                case 3: return SecondaryStructure::HELIX_PI;     // Pi helix
                case 5: return SecondaryStructure::HELIX_310;    // 3-10 helix
                default: return SecondaryStructure::HELIX_ALPHA;
            }
        }
    }
    
    // Check sheets
    for (const Sheet& sheet : cachedMoleculeData->sheets) {
        if (sheet.chain == chainId && 
            residueNum >= sheet.residueStart && 
            residueNum <= sheet.residueEnd) {
            return SecondaryStructure::SHEET;
        }
    }
    
    return SecondaryStructure::COIL;
}

void RepresentationRenderer::extractRibbonControlPoints(
    std::vector<RibbonControlPoint>& controlPoints
) {
    controlPoints.clear();

    if (!cachedMoleculeData) {
        std::cout << "[DEBUG] cachedMoleculeData is NULL!" << std::endl;
        return;
    }

    const auto& atoms = cachedMoleculeData->atoms;
    std::cout << "[DEBUG] Total atoms: " << atoms.size() << std::endl;

    // Print first 20 atom names RAW
    std::cout << "[DEBUG] First 20 atoms:" << std::endl;
    for (size_t i = 0; i < std::min(atoms.size(), (size_t)20); ++i) {
        std::cout << "  [" << i << "] name='" << atoms[i].name
            << "' res=" << atoms[i].residueName
            << " chain=" << atoms[i].chain
            << " resNum=" << atoms[i].residueNum << std::endl;
    }

    // Try to find ANY atom with "CA" anywhere in name
    int foundCA = 0;
    for (const auto& atom : atoms) {
        if (atom.name.find("CA") != std::string::npos) {
            foundCA++;
            if (foundCA <= 5) {
                std::cout << "[DEBUG] Found CA-like: '" << atom.name << "'" << std::endl;
            }
        }
    }
    std::cout << "[DEBUG] Atoms containing 'CA': " << foundCA << std::endl;

    // Build maps: (chain, residueNum) -> atom index for CA, CB, N, C atoms
    std::unordered_map<std::string, size_t> caAtomMap;
    std::unordered_map<std::string, size_t> cbAtomMap;
    std::unordered_map<std::string, size_t> nAtomMap;
    std::unordered_map<std::string, size_t> cAtomMap;

    for (size_t i = 0; i < atoms.size(); ++i) {
        const Atom& atom = atoms[i];
        std::string key = std::string(1, atom.chain) + "_" + std::to_string(atom.residueNum);

        // Get atom name and trim spaces
        std::string atomName = atom.name;
        while (!atomName.empty() && atomName[0] == ' ') atomName.erase(0, 1);
        while (!atomName.empty() && atomName.back() == ' ') atomName.pop_back();

        if (atomName == "CA") {
            caAtomMap[key] = i;
        }
        else if (atomName == "CB") {
            cbAtomMap[key] = i;
        }
        else if (atomName == "N") {
            nAtomMap[key] = i;
        }
        else if (atomName == "C") {
            cAtomMap[key] = i;
        }
    }

    std::cout << "[DEBUG] CA atoms mapped: " << caAtomMap.size() << std::endl;
    std::cout << "[DEBUG] CB atoms mapped: " << cbAtomMap.size() << std::endl;
    std::cout << "[DEBUG] N atoms mapped: " << nAtomMap.size() << std::endl;
    std::cout << "[DEBUG] C atoms mapped: " << cAtomMap.size() << std::endl;

    // Build residue list from sequence or atoms
    std::vector<std::pair<char, int>> residueList;

    // SKIP sequence - build directly from atoms instead
    // The sequence residue numbers often don't match PDB atom residue numbers
    std::cout << "[DEBUG] Building residue list from atoms..." << std::endl;
    std::set<std::pair<char, int>> seen;
    for (const Atom& atom : atoms) {
        auto key = std::make_pair(atom.chain, atom.residueNum);
        if (seen.find(key) == seen.end()) {
            seen.insert(key);
            residueList.push_back(key);
        }
    }
    std::sort(residueList.begin(), residueList.end());
    std::cout << "[DEBUG] Residues from atoms: " << residueList.size() << std::endl;

    // Now build control points
    char prevChain = '\0';
    int addedPoints = 0;

    for (const auto& [chainId, residueNum] : residueList) {
        std::string key = std::string(1, chainId) + "_" + std::to_string(residueNum);

        auto caIt = caAtomMap.find(key);
        if (caIt == caAtomMap.end()) continue;  // No CA atom for this residue

        const Atom& caAtom = atoms[caIt->second];

        RibbonControlPoint cp;
        cp.x = caAtom.coords.getX();
        cp.y = caAtom.coords.getY();
        cp.z = caAtom.coords.getZ();
        cp.chainId = chainId;
        cp.residueIndex = residueNum;
        cp.ssType = getSecondaryStructureForResidue(chainId, residueNum);

        // Compute guide normal from CA->CB direction
        auto cbIt = cbAtomMap.find(key);
        if (cbIt != cbAtomMap.end()) {
            const Atom& cbAtom = atoms[cbIt->second];
            cp.nx = cbAtom.coords.getX() - cp.x;
            cp.ny = cbAtom.coords.getY() - cp.y;
            cp.nz = cbAtom.coords.getZ() - cp.z;
        }
        else {
            // For glycine (no CB), use N->C perpendicular
            auto nIt = nAtomMap.find(key);
            auto cIt = cAtomMap.find(key);
            if (nIt != nAtomMap.end() && cIt != cAtomMap.end()) {
                const Atom& nAtom = atoms[nIt->second];
                const Atom& cAtom = atoms[cIt->second];
                float dx = cAtom.coords.getX() - nAtom.coords.getX();
                float dy = cAtom.coords.getY() - nAtom.coords.getY();
                float dz = cAtom.coords.getZ() - nAtom.coords.getZ();
                cp.nx = -dy;
                cp.ny = dx;
                cp.nz = 0.0f;
            }
            else {
                cp.nx = 0.0f;
                cp.ny = 1.0f;
                cp.nz = 0.0f;
            }
        }

        // Normalize
        float len = std::sqrt(cp.nx * cp.nx + cp.ny * cp.ny + cp.nz * cp.nz);
        if (len > 1e-6f) {
            cp.nx /= len;
            cp.ny /= len;
            cp.nz /= len;
        }

        // Set color based on secondary structure
        getSSColor(cp.ssType, cp.colorR, cp.colorG, cp.colorB);
        cp.colorA = 1.0f;

        prevChain = chainId;
        controlPoints.push_back(cp);
        addedPoints++;
    }

    std::cout << "[DEBUG] Final control points: " << addedPoints << std::endl;
}

void RepresentationRenderer::extractSurfaceAtoms(std::vector<SurfaceAtom>& surfaceAtoms) {
    surfaceAtoms.clear();
    if (!cachedMoleculeData) return;
    
    const auto& atoms = cachedMoleculeData->atoms;
    surfaceAtoms.reserve(atoms.size());
    
    for (const Atom& atom : atoms) {
        SurfaceAtom sa;
        sa.x = atom.coords.getX();
        sa.y = atom.coords.getY();
        sa.z = atom.coords.getZ();
        sa.radius = getVDWRadius(atom.element);
        getCPKColor(atom.element, sa.r, sa.g, sa.b);
        sa.a = 1.0f;
        sa.residueIndex = atom.residueNum;
        sa.chainId = atom.chain;
        
        surfaceAtoms.push_back(sa);
    }
}

void RepresentationRenderer::uploadRibbonBuffers() {
    if (ribbonTemplate.vertexCount() == 0) {
        ribbonBuffersValid = false;
        return;
    }
    
    glBindVertexArray(ribbonVAO);
    
    // Upload vertices
    glBindBuffer(GL_ARRAY_BUFFER, ribbonVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 ribbonTemplate.vertexCount() * RibbonTemplate::VERTEX_STRIDE,
                 ribbonTemplate.getVertexData(),
                 GL_STATIC_DRAW);
    
    // Upload indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ribbonEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 ribbonTemplate.indexCount() * sizeof(unsigned int),
                 ribbonTemplate.getIndexData(),
                 GL_STATIC_DRAW);
    
    // Set up vertex attributes
    // Position (location 0)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, RibbonTemplate::VERTEX_STRIDE,
                          (void*)offsetof(RibbonVertex, x));
    glEnableVertexAttribArray(0);
    
    // Normal (location 1)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, RibbonTemplate::VERTEX_STRIDE,
                          (void*)offsetof(RibbonVertex, nx));
    glEnableVertexAttribArray(1);
    
    // UV (location 2)
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, RibbonTemplate::VERTEX_STRIDE,
                          (void*)offsetof(RibbonVertex, u));
    glEnableVertexAttribArray(2);
    
    // Color (location 3)
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, RibbonTemplate::VERTEX_STRIDE,
                          (void*)offsetof(RibbonVertex, r));
    glEnableVertexAttribArray(3);
    
    glBindVertexArray(0);
    
    ribbonBuffersValid = true;
}

void RepresentationRenderer::uploadSurfaceBuffers() {
    if (surfaceTemplate.vertexCount() == 0) {
        surfaceBuffersValid = false;
        return;
    }
    
    glBindVertexArray(surfaceVAO);
    
    // Upload vertices
    glBindBuffer(GL_ARRAY_BUFFER, surfaceVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 surfaceTemplate.vertexCount() * SurfaceTemplate::VERTEX_STRIDE,
                 surfaceTemplate.getVertexData(),
                 GL_STATIC_DRAW);
    
    // Upload indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 surfaceTemplate.indexCount() * sizeof(unsigned int),
                 surfaceTemplate.getIndexData(),
                 GL_STATIC_DRAW);
    
    // Set up vertex attributes
    // Position (location 0)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, SurfaceTemplate::VERTEX_STRIDE,
                          (void*)offsetof(SurfaceVertex, x));
    glEnableVertexAttribArray(0);
    
    // Normal (location 1)
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, SurfaceTemplate::VERTEX_STRIDE,
                          (void*)offsetof(SurfaceVertex, nx));
    glEnableVertexAttribArray(1);
    
    // Color (location 2)
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, SurfaceTemplate::VERTEX_STRIDE,
                          (void*)offsetof(SurfaceVertex, r));
    glEnableVertexAttribArray(2);
    
    glBindVertexArray(0);
    
    surfaceBuffersValid = true;
}

void RepresentationRenderer::loadMoleculeData(MoleculeData* moleculeData) {
    std::cout << "[REPR] loadMoleculeData called" << std::endl;

    cachedMoleculeData = moleculeData;
    ribbonBuffersValid = false;
    surfaceBuffersValid = false;

    if (!moleculeData) {
        std::cout << "[REPR] moleculeData is null, returning" << std::endl;
        return;
    }

    std::cout << "[REPR] Atom count: " << moleculeData->atoms.size() << std::endl;
    std::cout << "[REPR] Starting ribbon generation..." << std::endl;

    regenerateRibbon();

    std::cout << "[REPR] Ribbon complete!" << std::endl;

    regenerateSurface();
}

void RepresentationRenderer::regenerateRibbon() {
    std::cout << "[REPR] extractRibbonControlPoints..." << std::endl;

    std::vector<RibbonControlPoint> controlPoints;
    extractRibbonControlPoints(controlPoints);

    std::cout << "[REPR] Control points: " << controlPoints.size() << std::endl;
    std::cout << "[REPR] generateRibbon..." << std::endl;

    ribbonTemplate.generateRibbon(controlPoints);

    std::cout << "[REPR] uploadRibbonBuffers..." << std::endl;

    uploadRibbonBuffers();

    std::cout << "[REPR] regenerateRibbon done!" << std::endl;
}

void RepresentationRenderer::regenerateSurface() {
    std::cout << "[REPR] start regenerateSurface..." << std::endl;

    std::vector<SurfaceAtom> surfaceAtoms;
    extractSurfaceAtoms(surfaceAtoms);
    surfaceTemplate.generateSurface(surfaceAtoms, surfaceType);
    uploadSurfaceBuffers();
    std::cout << "[REPR] regenerateSurface done!" << std::endl;

}

void RepresentationRenderer::setRepresentationEnabled(RepresentationType type, bool enabled) {
    RepresentationFlags flag = RepresentationFlags::NONE;
    
    switch (type) {
        case RepresentationType::BALL_AND_STICK:
            flag = RepresentationFlags::BALL_AND_STICK;
            break;
        case RepresentationType::RIBBON:
            flag = RepresentationFlags::RIBBON;
            break;
        case RepresentationType::SURFACE:
            flag = RepresentationFlags::SURFACE;
            break;
        case RepresentationType::SPACEFILL:
            flag = RepresentationFlags::SPACEFILL;
            break;
        case RepresentationType::WIREFRAME:
            flag = RepresentationFlags::WIREFRAME;
            break;
        case RepresentationType::STICK:
            flag = RepresentationFlags::STICK;
            break;
        default:
            return;
    }
    
    if (enabled) {
        activeRepresentations |= flag;
    } else {
        activeRepresentations = static_cast<RepresentationFlags>(
            static_cast<unsigned int>(activeRepresentations) & 
            ~static_cast<unsigned int>(flag)
        );
    }
}

bool RepresentationRenderer::isRepresentationEnabled(RepresentationType type) const {
    switch (type) {
        case RepresentationType::BALL_AND_STICK:
            return hasFlag(activeRepresentations, RepresentationFlags::BALL_AND_STICK);
        case RepresentationType::RIBBON:
            return hasFlag(activeRepresentations, RepresentationFlags::RIBBON);
        case RepresentationType::SURFACE:
            return hasFlag(activeRepresentations, RepresentationFlags::SURFACE);
        case RepresentationType::SPACEFILL:
            return hasFlag(activeRepresentations, RepresentationFlags::SPACEFILL);
        case RepresentationType::WIREFRAME:
            return hasFlag(activeRepresentations, RepresentationFlags::WIREFRAME);
        case RepresentationType::STICK:
            return hasFlag(activeRepresentations, RepresentationFlags::STICK);
        default:
            return false;
    }
}

void RepresentationRenderer::toggleRepresentation(RepresentationType type) {
    setRepresentationEnabled(type, !isRepresentationEnabled(type));
}

void RepresentationRenderer::setActiveRepresentations(RepresentationFlags flags) {
    activeRepresentations = flags;
}

RepresentationFlags RepresentationRenderer::getActiveRepresentations() const {
    return activeRepresentations;
}

void RepresentationRenderer::setSurfaceType(SurfaceTemplate::SurfaceType type) {
    if (surfaceType != type) {
        surfaceType = type;
        regenerateSurface();
    }
}

size_t RepresentationRenderer::getRibbonVertexCount() const {
    return ribbonTemplate.vertexCount();
}

size_t RepresentationRenderer::getSurfaceVertexCount() const {
    return surfaceTemplate.vertexCount();
}

void RepresentationRenderer::renderRibbon(
    Shader* shader,
    const Window* window,
    const Camera* camera,
    const float* modelMatrix
) {
    if (!ribbonBuffersValid || ribbonTemplate.indexCount() == 0) return;
    
    // TODO: Use your actual Shader interface
    // shader->use();
    // shader->setMat4("u_model", modelMatrix);
    // shader->setMat4("u_view", camera->getViewMatrix());
    // shader->setMat4("u_projection", camera->getProjectionMatrix(window));
    // shader->setVec3("u_lightDir", 0.3f, 0.5f, 0.8f);
    // shader->setVec3("u_cameraPos", camera->getPosition());
    
    glBindVertexArray(ribbonVAO);
    glDrawElements(GL_TRIANGLES, 
                   static_cast<GLsizei>(ribbonTemplate.indexCount()),
                   GL_UNSIGNED_INT, 
                   nullptr);
    glBindVertexArray(0);
}

void RepresentationRenderer::renderSurface(
    Shader* shader,
    const Window* window,
    const Camera* camera,
    const float* modelMatrix
) {
    if (!surfaceBuffersValid || surfaceTemplate.indexCount() == 0) return;
    
    // Enable blending for transparent surface
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // TODO: Use your actual Shader interface
    // shader->use();
    // shader->setMat4("u_model", modelMatrix);
    // shader->setMat4("u_view", camera->getViewMatrix());
    // shader->setMat4("u_projection", camera->getProjectionMatrix(window));
    // shader->setVec3("u_lightDir", 0.3f, 0.5f, 0.8f);
    // shader->setVec3("u_cameraPos", camera->getPosition());
    // shader->setFloat("u_opacity", surfaceOpacity);
    
    glBindVertexArray(surfaceVAO);
    glDrawElements(GL_TRIANGLES,
                   static_cast<GLsizei>(surfaceTemplate.indexCount()),
                   GL_UNSIGNED_INT,
                   nullptr);
    glBindVertexArray(0);
    
    glDisable(GL_BLEND);
}

void RepresentationRenderer::render(
    Shader* ribbonShader,
    Shader* surfaceShader,
    const Window* window,
    const Camera* camera,
    const float* modelMatrix
) {
    // Ball-and-stick is rendered by existing Model class
    // This renderer handles ribbon and surface
    
    if (hasFlag(activeRepresentations, RepresentationFlags::RIBBON)) {
        renderRibbon(ribbonShader, window, camera, modelMatrix);
    }
    
    if (hasFlag(activeRepresentations, RepresentationFlags::SURFACE)) {
        renderSurface(surfaceShader, window, camera, modelMatrix);
    }
}
