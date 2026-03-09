#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <optional>
#include <fstream>
#include <iostream>
#include "Parser.h"

// Represents one alignment segment between PDB and UniProt
struct AlignmentSegment {
    int pdbStart;      // PDB_BEG
    int pdbEnd;        // PDB_END
    int uniprotStart;  // SP_BEG
    int uniprotEnd;    // SP_END
    std::string uniprotID;  // SP_PRIMARY

    // Check if a PDB position falls within this segment
    bool contains(int pdbPosition) const {
        return pdbPosition >= pdbStart && pdbPosition <= pdbEnd;
    }

    // Convert PDB position to UniProt position
    int toUniProtPosition(int pdbPosition) const {
        int offset = pdbPosition - pdbStart;
        return uniprotStart + offset;
    }
};

// Stores all alignment info for one PDB chain
struct ChainMapping {
    std::string pdbID;
    char chain;
    std::string uniprotID;
    std::vector<AlignmentSegment> segments;
};

class PDBUniProtMapper {
private:
    // Key: "PDBID_CHAIN" (e.g., "101m_A") -> ChainMapping
    std::unordered_map<std::string, ChainMapping> mappings;
    bool isLoaded = false;

    // Create lookup key
    std::string makeKey(const std::string& pdbID, char chain) const {
        return pdbID + "_" + std::string(1, chain);
    }

public:
    // Load mappings from CSV file
    bool loadFromCSV(const std::string& filepath);

    // Convert PDB position to UniProt position
    std::optional<int> pdbToUniProtPosition(
        const std::string& pdbID,
        char chain,
        int pdbPosition
    ) const;

    // Get UniProt ID for a PDB chain
    std::optional<std::string> getUniProtID(
        const std::string& pdbID,
        char chain
    ) const;

    // Get full mapping info for a chain
    std::optional<ChainMapping> getChainMapping(
        const std::string& pdbID,
        char chain
    ) const;

    bool isReady() const { return isLoaded; }
    size_t size() const { return mappings.size(); }
    std::vector<ChainMapping> getAllChainsForPDB(const std::string& pdbID) const;
    void displayPDBInfo(const std::string& pdbID) const;

};

struct PDBUniProtMapping {
    std::string pdbID;
    char chain;
    std::string uniprotID;
    int pdbStart;
    int pdbEnd;
    int uniprotStart;
    int uniprotEnd;

    // Add constructor
    PDBUniProtMapping() = default;

    PDBUniProtMapping(
        std::string pdb,
        char ch,
        std::string uniprot,
        int pStart,
        int pEnd,
        int uStart,
        int uEnd
    ) : pdbID(std::move(pdb)),
        chain(ch),
        uniprotID(std::move(uniprot)),
        pdbStart(pStart),
        pdbEnd(pEnd),
        uniprotStart(uStart),
        uniprotEnd(uEnd) {
    }
};