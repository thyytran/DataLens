#pragma once
#include <alphamissense/AlphaMissenseDB.h>
#include <bio/MoleculeData.h>
#include <mapper/PDBToUniProtMapper.h>

struct MutationResult {
    // Input
    char chain;
    int position;
    char originalAA;
    char mutantAA;
    std::string uniprotID;

    // AlphaMissense results
    float pathogenicityScore;
    std::string pathogenicityClass;

    void print() const;
};

class MutationAnalyzer {
private:
    AlphaMissenseDB* alphaMissenseDB;
    PDBUniProtMapper* pdbMapper;  // Add this

    // Get UniProt ID
    std::string getUniProtID(const std::string& pdbID, char chain);

    // Extract original AA from structure
    char getOriginalAA(MoleculeData* moleculeData, char chain, int position);

public:
    MutationAnalyzer(AlphaMissenseDB* db, PDBUniProtMapper* mapper) : alphaMissenseDB(db), pdbMapper(mapper) {}

    // Main analysis function
    std::optional<MutationResult> analyzeMutation(
        const std::string& pdbID,
        MoleculeData* moleculeData,
        char chain,
        int position,
        char mutantAA
    );
};