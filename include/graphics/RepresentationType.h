#pragma once

/*
 * RepresentationType.h
 * Defines the different molecular representation types supported by DataLens.
 * 
 * Three primary representations:
 * 1. Ball-and-Stick: Atoms as spheres, bonds as cylinders (existing)
 * 2. Ribbon/Cartoon: Secondary structure visualization
 * 3. Surface: Molecular surface (VDW, SAS, SES)
 */

enum class RepresentationType {
    NONE = 0,
    BALL_AND_STICK = 1,    // Atoms as spheres, bonds as cylinders
    RIBBON = 2,            // Cartoon/ribbon for secondary structure
    SURFACE = 3,           // Molecular surface
    SPACEFILL = 4,         // VDW spheres (full radius, no bonds)
    WIREFRAME = 5,         // Lines only
    STICK = 6              // Bonds only, no atom spheres
};

// Bit flags for combining multiple representations
enum class RepresentationFlags : unsigned int {
    NONE           = 0,
    BALL_AND_STICK = 1 << 0,
    RIBBON         = 1 << 1,
    SURFACE        = 1 << 2,
    SPACEFILL      = 1 << 3,
    WIREFRAME      = 1 << 4,
    STICK          = 1 << 5
};

// Bitwise operators for RepresentationFlags
inline RepresentationFlags operator|(RepresentationFlags a, RepresentationFlags b) {
    return static_cast<RepresentationFlags>(
        static_cast<unsigned int>(a) | static_cast<unsigned int>(b)
    );
}

inline RepresentationFlags operator&(RepresentationFlags a, RepresentationFlags b) {
    return static_cast<RepresentationFlags>(
        static_cast<unsigned int>(a) & static_cast<unsigned int>(b)
    );
}

inline RepresentationFlags& operator|=(RepresentationFlags& a, RepresentationFlags b) {
    a = a | b;
    return a;
}

inline RepresentationFlags& operator&=(RepresentationFlags& a, RepresentationFlags b) {
    a = a & b;
    return a;
}

inline bool hasFlag(RepresentationFlags flags, RepresentationFlags flag) {
    return (static_cast<unsigned int>(flags) & static_cast<unsigned int>(flag)) != 0;
}
