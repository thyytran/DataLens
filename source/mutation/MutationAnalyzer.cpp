// MutationAnalyzer.cpp
#include "mutation/MutationAnalyzer.h"
#include "HTTPConnection.h"
#include "Parser.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <unordered_map>

using json = nlohmann::json;

// Helper: Convert 3-letter to 1-letter
namespace {
    const std::unordered_map<std::string, char> threeToOne = {
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'},
        {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'},
        {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'},
        {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'}, {"SER", 'S'},
        {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'}
    };
}

std::string MutationAnalyzer::getUniProtID(const std::string& pdbID, char chain) {
    HTTPConnection conn;
    std::string url = "https://data.rcsb.org/rest/v1/core/uniprot/" +
        pdbID + "/" + std::string(1, chain);

    if (conn.get(url)) {
        try {
            json response = json::parse(conn.response);
            return response["rcsb_uniprot_container_identifiers"][0]
                ["uniprot_id"].get<std::string>();
        }
        catch (...) {
            std::cerr << "  Error parsing UniProt response\n";
        }
    }
    return "";
}

char MutationAnalyzer::getOriginalAA(MoleculeData* moleculeData, char chain, int position) {
    if (!moleculeData) return 'X';

    for (const auto& atom : moleculeData->atoms) {
        if (atom.chain == chain && atom.residueNum == position) {
            auto it = threeToOne.find(atom.residueName);
            if (it != threeToOne.end()) {
                return it->second;
            }
        }
    }
    return 'X';
}

void MutationResult::print() const {
    std::cout << "\n================================================\n";
    std::cout << "        Mutation Analysis Results              \n";
    std::cout << "================================================\n\n";

    std::cout << "Mutation: Chain " << chain << " - "
        << originalAA << position << mutantAA << "\n";
    std::cout << "UniProt ID: " << uniprotID << "\n";
    std::cout << "------------------------------------------------\n\n";

    std::cout << "AlphaMissense Pathogenicity Score\n";
    if (pathogenicityScore >= 0) {
        std::cout << "   Score: " << pathogenicityScore << "\n";
        std::cout << "   Classification: " << pathogenicityClass << "\n\n";

        std::cout << "Interpretation:\n";
        if (pathogenicityScore > 0.8) {
            std::cout << "   >> High confidence PATHOGENIC mutation\n";
            std::cout << "   This variant is very likely to cause disease.\n";
        }
        else if (pathogenicityScore > 0.564) {
            std::cout << "   >> Likely PATHOGENIC mutation\n";
            std::cout << "   This variant is predicted to be disease-causing.\n";
        }
        else if (pathogenicityScore > 0.34) {
            std::cout << "   >> UNCERTAIN significance\n";
            std::cout << "   Effect is ambiguous - requires further study.\n";
        }
        else {
            std::cout << "   >> Likely BENIGN mutation\n";
            std::cout << "   This variant is predicted to be tolerated.\n";
        }
    }
    else {
        std::cout << "   Not available in AlphaMissense database\n";
        std::cout << "   This protein/position may not be covered.\n";
    }

    std::cout << "\n================================================\n\n";
}

std::optional<MutationResult> MutationAnalyzer::analyzeMutation(
    const std::string& pdbID,
    MoleculeData* moleculeData,
    char chain,
    int pdbPosition,
    char mutantAA
) {
    std::cout << "\nAnalyzing Mutation...\n";
    std::cout << "   PDB: " << pdbID << "\n";
    std::cout << "   Chain: " << chain << "\n";
    std::cout << "   PDB Position: " << pdbPosition << "\n";
    std::cout << "   Target AA: " << mutantAA << "\n\n";

    MutationResult result;
    result.chain = chain;
    result.position = pdbPosition;
    result.mutantAA = mutantAA;

    // Get original amino acid
    std::cout << "   Extracting original amino acid...\n";
    result.originalAA = getOriginalAA(moleculeData, chain, pdbPosition);
    if (result.originalAA == 'X') {
        std::cerr << "Error: Could not find residue at position "
            << pdbPosition << " chain " << chain << "\n\n";
        return std::nullopt;
    }

    std::cout << "   Original AA: " << result.originalAA << "\n";

    // Check if mutation is to same amino acid
    if (result.originalAA == mutantAA) {
        std::cout << "Warning: Target is same as original. No mutation.\n\n";
        return std::nullopt;
    }

    // Get UniProt ID using mapper
    std::cout << "   Looking up UniProt ID...\n";
    auto uniprotID = pdbMapper->getUniProtID(pdbID, chain);
    if (!uniprotID) {
        std::cerr << "Error: No UniProt mapping found for "
            << pdbID << " chain " << chain << "\n\n";
        return std::nullopt;
    }
    result.uniprotID = *uniprotID;
    std::cout << "   UniProt ID: " << result.uniprotID << "\n";

    // Convert PDB position to UniProt position
    std::cout << "   Converting PDB position to UniProt position...\n";
    auto uniprotPosition = pdbMapper->pdbToUniProtPosition(pdbID, chain, pdbPosition);
    if (!uniprotPosition) {
        std::cerr << "Error: PDB position " << pdbPosition
            << " not found in UniProt mapping\n\n";
        return std::nullopt;
    }

    std::cout << "   UniProt Position: " << *uniprotPosition << "\n";
    std::cout << "   Mutation: " << result.originalAA << *uniprotPosition
        << mutantAA << "\n\n";

    // Query AlphaMissense
    std::cout << "   Querying AlphaMissense...\n";
    auto score = alphaMissenseDB->getScore(
        result.uniprotID,
        result.originalAA,
        *uniprotPosition,
        mutantAA
    );

    if (score) {
        result.pathogenicityScore = score->pathogenicity;
        result.pathogenicityClass = score->classification;
        std::cout << "   Found in AlphaMissense\n";
    }
    else {
        std::cout << "   Not found in AlphaMissense database\n";
        result.pathogenicityScore = -1.0f;
        result.pathogenicityClass = "unknown";
    }

    return result;
}