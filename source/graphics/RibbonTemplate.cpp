#include "RibbonTemplate.h"
#include <cmath>
#include <algorithm>

/*
 * RibbonTemplate.cpp
 * Implementation of ribbon/cartoon geometry generation.
 */

void RibbonTemplate::clear() {
    vertices.clear();
    indices.clear();
}

void RibbonTemplate::normalize(float& x, float& y, float& z) {
    float len = std::sqrt(x*x + y*y + z*z);
    if (len > 1e-6f) {
        x /= len;
        y /= len;
        z /= len;
    }
}

void RibbonTemplate::cross(
    float ax, float ay, float az,
    float bx, float by, float bz,
    float& rx, float& ry, float& rz
) {
    rx = ay * bz - az * by;
    ry = az * bx - ax * bz;
    rz = ax * by - ay * bx;
}

void RibbonTemplate::catmullRom(
    float p0x, float p0y, float p0z,
    float p1x, float p1y, float p1z,
    float p2x, float p2y, float p2z,
    float p3x, float p3y, float p3z,
    float t,
    float& outX, float& outY, float& outZ
) {
    float t2 = t * t;
    float t3 = t2 * t;
    
    // Catmull-Rom basis functions
    float b0 = -0.5f*t3 + t2 - 0.5f*t;
    float b1 = 1.5f*t3 - 2.5f*t2 + 1.0f;
    float b2 = -1.5f*t3 + 2.0f*t2 + 0.5f*t;
    float b3 = 0.5f*t3 - 0.5f*t2;
    
    outX = b0*p0x + b1*p1x + b2*p2x + b3*p3x;
    outY = b0*p0y + b1*p1y + b2*p2y + b3*p3y;
    outZ = b0*p0z + b1*p1z + b2*p2z + b3*p3z;
}

void RibbonTemplate::catmullRomTangent(
    float p0x, float p0y, float p0z,
    float p1x, float p1y, float p1z,
    float p2x, float p2y, float p2z,
    float p3x, float p3y, float p3z,
    float t,
    float& outX, float& outY, float& outZ
) {
    float t2 = t * t;
    
    // Derivatives of Catmull-Rom basis functions
    float b0 = -1.5f*t2 + 2.0f*t - 0.5f;
    float b1 = 4.5f*t2 - 5.0f*t;
    float b2 = -4.5f*t2 + 4.0f*t + 0.5f;
    float b3 = 1.5f*t2 - t;
    
    outX = b0*p0x + b1*p1x + b2*p2x + b3*p3x;
    outY = b0*p0y + b1*p1y + b2*p2y + b3*p3y;
    outZ = b0*p0z + b1*p1z + b2*p2z + b3*p3z;
    
    normalize(outX, outY, outZ);
}

void RibbonTemplate::generateHelixCrossSection(
    float cx, float cy, float cz,
    float tx, float ty, float tz,
    float nx, float ny, float nz,
    float width, float thickness,
    float r, float g, float b, float a,
    float v,
    std::vector<RibbonVertex>& crossSection
) {
    // Compute binormal (perpendicular to tangent and normal)
    float bx, by, bz;
    cross(tx, ty, tz, nx, ny, nz, bx, by, bz);
    normalize(bx, by, bz);
    
    // Recompute normal to ensure orthogonality
    cross(bx, by, bz, tx, ty, tz, nx, ny, nz);
    normalize(nx, ny, nz);
    
    crossSection.clear();
    
    // 4-point rectangular cross section (flat ribbon)
    // Top-left, top-right, bottom-right, bottom-left
    float hw = width * 0.5f;
    float ht = thickness * 0.5f;
    
    // Top-left
    RibbonVertex v0;
    v0.x = cx - hw * bx + ht * nx;
    v0.y = cy - hw * by + ht * ny;
    v0.z = cz - hw * bz + ht * nz;
    v0.nx = nx; v0.ny = ny; v0.nz = nz;  // Top face normal
    v0.u = 0.0f; v0.v = v;
    v0.r = r; v0.g = g; v0.b = b; v0.a = a;
    crossSection.push_back(v0);
    
    // Top-right
    RibbonVertex v1;
    v1.x = cx + hw * bx + ht * nx;
    v1.y = cy + hw * by + ht * ny;
    v1.z = cz + hw * bz + ht * nz;
    v1.nx = nx; v1.ny = ny; v1.nz = nz;
    v1.u = 1.0f; v1.v = v;
    v1.r = r; v1.g = g; v1.b = b; v1.a = a;
    crossSection.push_back(v1);
    
    // Bottom-right
    RibbonVertex v2;
    v2.x = cx + hw * bx - ht * nx;
    v2.y = cy + hw * by - ht * ny;
    v2.z = cz + hw * bz - ht * nz;
    v2.nx = -nx; v2.ny = -ny; v2.nz = -nz;  // Bottom face normal
    v2.u = 1.0f; v2.v = v;
    v2.r = r; v2.g = g; v2.b = b; v2.a = a;
    crossSection.push_back(v2);
    
    // Bottom-left
    RibbonVertex v3;
    v3.x = cx - hw * bx - ht * nx;
    v3.y = cy - hw * by - ht * ny;
    v3.z = cz - hw * bz - ht * nz;
    v3.nx = -nx; v3.ny = -ny; v3.nz = -nz;
    v3.u = 0.0f; v3.v = v;
    v3.r = r; v3.g = g; v3.b = b; v3.a = a;
    crossSection.push_back(v3);
}

void RibbonTemplate::generateSheetCrossSection(
    float cx, float cy, float cz,
    float tx, float ty, float tz,
    float nx, float ny, float nz,
    float width, float thickness,
    float r, float g, float b, float a,
    float v,
    bool isArrowHead,
    std::vector<RibbonVertex>& crossSection
) {
    // Use wider width for arrow head
    float actualWidth = isArrowHead ? SHEET_ARROW_WIDTH : width;
    generateHelixCrossSection(cx, cy, cz, tx, ty, tz, nx, ny, nz,
                              actualWidth, thickness, r, g, b, a, v, crossSection);
}

void RibbonTemplate::generateCoilCrossSection(
    float cx, float cy, float cz,
    float tx, float ty, float tz,
    float nx, float ny, float nz,
    float radius,
    float r, float g, float b, float a,
    float v,
    std::vector<RibbonVertex>& crossSection
) {
    // Compute binormal
    float bx, by, bz;
    cross(tx, ty, tz, nx, ny, nz, bx, by, bz);
    normalize(bx, by, bz);
    
    // Recompute normal
    cross(bx, by, bz, tx, ty, tz, nx, ny, nz);
    normalize(nx, ny, nz);
    
    crossSection.clear();
    
    // Circular cross section
    float angleStep = 2.0f * 3.14159265f / TUBE_SIDES;
    for (int i = 0; i < TUBE_SIDES; ++i) {
        float angle = i * angleStep;
        float cosA = std::cos(angle);
        float sinA = std::sin(angle);
        
        // Direction from center
        float dirX = cosA * nx + sinA * bx;
        float dirY = cosA * ny + sinA * by;
        float dirZ = cosA * nz + sinA * bz;
        
        RibbonVertex vert;
        vert.x = cx + radius * dirX;
        vert.y = cy + radius * dirY;
        vert.z = cz + radius * dirZ;
        vert.nx = dirX;
        vert.ny = dirY;
        vert.nz = dirZ;
        vert.u = static_cast<float>(i) / TUBE_SIDES;
        vert.v = v;
        vert.r = r; vert.g = g; vert.b = b; vert.a = a;
        
        crossSection.push_back(vert);
    }
}

void RibbonTemplate::connectCrossSections(
    const std::vector<RibbonVertex>& cs1,
    const std::vector<RibbonVertex>& cs2,
    unsigned int baseIndex1,
    unsigned int baseIndex2
) {
    size_t n = cs1.size();
    if (n != cs2.size() || n < 2) return;
    
    for (size_t i = 0; i < n; ++i) {
        size_t next = (i + 1) % n;
        
        unsigned int i1 = baseIndex1 + static_cast<unsigned int>(i);
        unsigned int i2 = baseIndex1 + static_cast<unsigned int>(next);
        unsigned int i3 = baseIndex2 + static_cast<unsigned int>(i);
        unsigned int i4 = baseIndex2 + static_cast<unsigned int>(next);
        
        // Two triangles per quad
        // Triangle 1: i1, i3, i2
        indices.push_back(i1);
        indices.push_back(i3);
        indices.push_back(i2);
        
        // Triangle 2: i2, i3, i4
        indices.push_back(i2);
        indices.push_back(i3);
        indices.push_back(i4);
    }
}

void RibbonTemplate::generateRibbon(const std::vector<RibbonControlPoint>& controlPoints) {
    clear();
    
    if (controlPoints.size() < 2) return;
    
    std::vector<RibbonVertex> prevCrossSection;
    unsigned int prevBaseIndex = 0;
    
    int numPoints = static_cast<int>(controlPoints.size());
    
    for (int i = 0; i < numPoints - 1; ++i) {
        // Get four points for Catmull-Rom (clamp at boundaries)
        int i0 = std::max(0, i - 1);
        int i1 = i;
        int i2 = i + 1;
        int i3 = std::min(numPoints - 1, i + 2);
        
        const auto& p0 = controlPoints[i0];
        const auto& p1 = controlPoints[i1];
        const auto& p2 = controlPoints[i2];
        const auto& p3 = controlPoints[i3];
        
        // Check for chain break
        if (p1.chainId != p2.chainId) {
            prevCrossSection.clear();
            continue;
        }
        
        // Generate segments between p1 and p2
        for (int s = 0; s <= SPLINE_SEGMENTS; ++s) {
            float t = static_cast<float>(s) / SPLINE_SEGMENTS;
            
            // Skip first segment if we already have a previous cross-section
            if (s == 0 && !prevCrossSection.empty()) continue;
            
            // Interpolate position
            float px, py, pz;
            catmullRom(p0.x, p0.y, p0.z,
                       p1.x, p1.y, p1.z,
                       p2.x, p2.y, p2.z,
                       p3.x, p3.y, p3.z,
                       t, px, py, pz);
            
            // Interpolate tangent
            float tx, ty, tz;
            catmullRomTangent(p0.x, p0.y, p0.z,
                              p1.x, p1.y, p1.z,
                              p2.x, p2.y, p2.z,
                              p3.x, p3.y, p3.z,
                              t, tx, ty, tz);
            
            // Interpolate normal
            float nx = (1.0f - t) * p1.nx + t * p2.nx;
            float ny = (1.0f - t) * p1.ny + t * p2.ny;
            float nz = (1.0f - t) * p1.nz + t * p2.nz;
            normalize(nx, ny, nz);
            
            // Interpolate color
            float r = (1.0f - t) * p1.colorR + t * p2.colorR;
            float g = (1.0f - t) * p1.colorG + t * p2.colorG;
            float b = (1.0f - t) * p1.colorB + t * p2.colorB;
            float a = (1.0f - t) * p1.colorA + t * p2.colorA;
            
            // Determine secondary structure type (use p1's type for first half, p2's for second)
            SecondaryStructure ssType = (t < 0.5f) ? p1.ssType : p2.ssType;
            
            // V coordinate for texture mapping
            float v = static_cast<float>(i) + t;
            
            // Generate cross section based on secondary structure type
            std::vector<RibbonVertex> currentCrossSection;
            
            if (isHelix(ssType)) {
                generateHelixCrossSection(px, py, pz, tx, ty, tz, nx, ny, nz,
                                         HELIX_WIDTH, HELIX_THICKNESS,
                                         r, g, b, a, v, currentCrossSection);
            }
            else if (isSheet(ssType)) {
                // Check if near end of strand for arrow head
                bool isArrowHead = (i == numPoints - 2 && s >= SPLINE_SEGMENTS - 1) ||
                                   (p2.ssType != SecondaryStructure::SHEET && t > 0.7f);
                generateSheetCrossSection(px, py, pz, tx, ty, tz, nx, ny, nz,
                                         SHEET_WIDTH, SHEET_THICKNESS,
                                         r, g, b, a, v, isArrowHead, currentCrossSection);
            }
            else {
                generateCoilCrossSection(px, py, pz, tx, ty, tz, nx, ny, nz,
                                        COIL_RADIUS, r, g, b, a, v, currentCrossSection);
            }
            
            // Add vertices to main array
            unsigned int currentBaseIndex = static_cast<unsigned int>(vertices.size());
            for (const auto& vert : currentCrossSection) {
                vertices.push_back(vert);
            }
            
            // Connect to previous cross-section if compatible
            if (!prevCrossSection.empty() && prevCrossSection.size() == currentCrossSection.size()) {
                connectCrossSections(prevCrossSection, currentCrossSection,
                                    prevBaseIndex, currentBaseIndex);
            }
            
            prevCrossSection = currentCrossSection;
            prevBaseIndex = currentBaseIndex;
        }
    }
}
