#pragma once
#include <math/Vec.h>
#include <vector>

// Per-residue data (built once on PDB load)
enum class SSType { COIL, HELIX, SHEET };

struct ResidueInfo {
    Vec3 caPosition;      // alpha carbon position
    Vec3 oPosition;       // carbonyl oxygen (defines ribbon normal)
    SSType ssType;
    char chain;
    int seqNum;
};

// Per-chain spline + generated mesh
struct RibbonChain {
    std::vector<Vec3> splinePoints;    // Catmull-Rom through CAs
    std::vector<Vec3> normals;         // from CA -> O direction, smoothed
    std::vector<Vec3> tangents;        // spline derivatives
    std::vector<SSType> ssPerPoint;    // interpolated SS assignment

    // Generated mesh (uploaded to GPU)
    std::vector<float> vertices;       // pos + normal per ribbon vertex
    std::vector<unsigned int> indices;
};