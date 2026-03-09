#pragma once

#include <string>
#include <vector>
#include <unordered_set>
#include "bio/PDBFile.h"
#include <alphamissense/AlphaMissenseDB.h>

struct MutationSuggestion {
    char chain;
    int position;
    char originalAA;
    char suggestedAA;
    float pathogenicityScore;
    std::string classification;

    std::string toString() const {
        return "Chain " + std::string(1, chain) + " - " +
            std::string(1, originalAA) + std::to_string(position) +
            std::string(1, suggestedAA) +
            " (score: " + std::to_string(pathogenicityScore).substr(0, 4) +
            ", " + classification + ")";
    }
};

struct ChainInfo {
    char chainID;
    std::string uniprotID;
    std::string sequence;
    std::vector<int> positions;
};

class MutationSuggester {
private:
    AlphaMissenseDB* alphaMissenseDB;

    // Get UniProt ID from PDB
    std::string getUniProtID(const std::string& pdbID, char chain);

    // Extract all chains from molecule data
    std::vector<ChainInfo> extractChainInfo(
        const std::string& pdbID,
        MoleculeData* moleculeData
    );

    // Get ALL mutations for a UniProt ID from AlphaMissense
    std::vector<MutationSuggestion> getAllMutationsForUniProt(
        const std::string& uniprotID,
        char chain,
        const std::string& sequence,
        const std::vector<int>& positions,
        float minPathogenicity = 0.7
    );

public:
    MutationSuggester(AlphaMissenseDB* db) : alphaMissenseDB(db) {}

    // Main suggestion method - now chain-agnostic
    std::vector<MutationSuggestion> suggestMutations(
        const std::string& pdbID,
        MoleculeData* moleculeData,
        int maxSuggestions = 20
    );
};