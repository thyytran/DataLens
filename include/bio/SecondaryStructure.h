#pragma once

#include <string>

/*
 * SecondaryStructure.h
 * Defines protein secondary structure types used for ribbon/cartoon rendering.
 * Based on DSSP classification and PDB HELIX/SHEET records.
 */

enum class SecondaryStructure {
    COIL = 0,          // Random coil / loop
    HELIX_ALPHA = 1,   // Alpha helix (most common)
    HELIX_310 = 2,     // 3-10 helix (tighter)
    HELIX_PI = 3,      // Pi helix (wider)
    SHEET = 4,         // Beta sheet / strand
    TURN = 5,          // Beta turn
    BRIDGE = 6         // Beta bridge (isolated)
};

/*
 * Convert PDB HELIX record type to SecondaryStructure.
 * PDB helix types: 1=right-handed alpha, 3=pi, 5=3-10
 */
inline SecondaryStructure helixTypeFromPDB(int pdbHelixType) {
    switch (pdbHelixType) {
        case 1: return SecondaryStructure::HELIX_ALPHA;
        case 3: return SecondaryStructure::HELIX_PI;
        case 5: return SecondaryStructure::HELIX_310;
        default: return SecondaryStructure::HELIX_ALPHA;
    }
}

/*
 * Check if a secondary structure type is a helix variant.
 */
inline bool isHelix(SecondaryStructure ss) {
    return ss == SecondaryStructure::HELIX_ALPHA ||
           ss == SecondaryStructure::HELIX_310 ||
           ss == SecondaryStructure::HELIX_PI;
}

/*
 * Check if a secondary structure type is a sheet/strand.
 */
inline bool isSheet(SecondaryStructure ss) {
    return ss == SecondaryStructure::SHEET ||
           ss == SecondaryStructure::BRIDGE;
}

/*
 * Check if a secondary structure type is coil/loop.
 */
inline bool isCoil(SecondaryStructure ss) {
    return ss == SecondaryStructure::COIL ||
           ss == SecondaryStructure::TURN;
}

/*
 * Get human-readable name for secondary structure type.
 */
inline std::string getSecondaryStructureName(SecondaryStructure ss) {
    switch (ss) {
        case SecondaryStructure::COIL:        return "Coil";
        case SecondaryStructure::HELIX_ALPHA: return "Alpha Helix";
        case SecondaryStructure::HELIX_310:   return "3-10 Helix";
        case SecondaryStructure::HELIX_PI:    return "Pi Helix";
        case SecondaryStructure::SHEET:       return "Beta Sheet";
        case SecondaryStructure::TURN:        return "Turn";
        case SecondaryStructure::BRIDGE:      return "Beta Bridge";
        default:                              return "Unknown";
    }
}
