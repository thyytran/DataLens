// MutationSuggester.cpp
#include "mutation/MutationSuggester.h"
#include "HTTPConnection.h"
#include "nlohmann/json.hpp"
#include "Parser.h"
#include <algorithm>
#include <unordered_map>
#include <set>

using json = nlohmann::json;

// Helper: Convert 3-letter amino acid code to 1-letter code
namespace {
    const std::unordered_map<std::string, char> threeToOne = {
        {"ALA", 'A'}, {"ARG", 'R'}, {"ASN", 'N'}, {"ASP", 'D'},
        {"CYS", 'C'}, {"GLN", 'Q'}, {"GLU", 'E'}, {"GLY", 'G'},
        {"HIS", 'H'}, {"ILE", 'I'}, {"LEU", 'L'}, {"LYS", 'K'},
        {"MET", 'M'}, {"PHE", 'F'}, {"PRO", 'P'}, {"SER", 'S'},
        {"THR", 'T'}, {"TRP", 'W'}, {"TYR", 'Y'}, {"VAL", 'V'}
    };

    char convertThreeToOne(const std::string& threeLetter) {
        std::string upper = Parser::uppercase(threeLetter);
        auto it = threeToOne.find(upper);
        if (it != threeToOne.end()) {
            return it->second;
        }
        return 'X';  // Unknown amino acid
    }
}

std::string MutationSuggester::getUniProtID(const std::string& pdbID, char chain) {
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
            // Silent fail - this chain might not map to UniProt
        }
    }
    return "";
}

std::vector<ChainInfo> MutationSuggester::extractChainInfo(
    const std::string& pdbID,
    MoleculeData* moleculeData
) {
    std::vector<ChainInfo> chainInfos;

    if (!moleculeData) return chainInfos;

    // Find all unique chains
    std::set<char> uniqueChains;
    for (const auto& atom : moleculeData->atoms) {
        uniqueChains.insert(atom.chain);
    }

    std::cout << "Found " << uniqueChains.size() << " chain(s): ";
    for (char c : uniqueChains) {
        std::cout << c << " ";
    }
    std::cout << "\n";

    // Extract sequence and UniProt ID for each chain
    for (char chain : uniqueChains) {
        ChainInfo info;
        info.chainID = chain;

        std::cout << "\nProcessing chain " << chain << "...\n";

        // Build sequence from atoms for this chain
        int lastResNum = -1;
        for (const auto& atom : moleculeData->atoms) {
            if (atom.chain == chain && atom.residueNum != lastResNum) {
                char residueCode = convertThreeToOne(atom.residueName);
                if (residueCode != 'X') {
                    info.sequence += residueCode;
                    info.positions.push_back(atom.residueNum);
                    lastResNum = atom.residueNum;
                }
            }
        }

        if (info.sequence.empty()) {
            std::cout << "  Chain " << chain << " has no valid sequence, skipping.\n";
            continue;
        }

        std::cout << "  Sequence length: " << info.sequence.length() << " residues\n";
        std::cout << "  Fetching UniProt ID...\n";

        // Get UniProt ID
        info.uniprotID = getUniProtID(pdbID, chain);

        if (info.uniprotID.empty()) {
            std::cout << "  Could not retrieve UniProt ID for chain " << chain << "\n";
            continue;
        }

        std::cout << "  UniProt ID: " << info.uniprotID << "\n";
        chainInfos.push_back(info);
    }

    return chainInfos;
}

std::vector<MutationSuggestion> MutationSuggester::getAllMutationsForUniProt(
    const std::string& uniprotID,
    char chain,
    const std::string& sequence,
    const std::vector<int>& positions,
    float minPathogenicity
) {
    std::vector<MutationSuggestion> mutations;

    std::cout << "  Searching AlphaMissense for pathogenic mutations...\n";

    // For each position in the sequence
    for (size_t i = 0; i < sequence.length(); ++i) {
        char originalAA = sequence[i];
        int position = positions[i];

        // Try all 20 amino acids
        for (char altAA = 'A'; altAA <= 'Z'; ++altAA) {
            // Skip invalid amino acids and same as original
            if (altAA == originalAA) continue;
            if (convertThreeToOne(std::string(1, altAA)) == 'X') continue;

            // Query AlphaMissense
            auto score = alphaMissenseDB->getScore(uniprotID, originalAA, position, altAA);

            if (score && score->pathogenicity >= minPathogenicity) {
                mutations.push_back({
                    chain,
                    position,
                    originalAA,
                    altAA,
                    score->pathogenicity,
                    score->classification
                    });
            }
        }
    }

    std::cout << "  Found " << mutations.size() << " mutations with pathogenicity >= "
        << minPathogenicity << "\n";

    return mutations;
}

std::vector<MutationSuggestion> MutationSuggester::suggestMutations(
    const std::string& pdbID,
    MoleculeData* moleculeData,
    int maxSuggestions
) {
    std::vector<MutationSuggestion> allSuggestions;

    std::cout << "\n=== Analyzing protein for mutation suggestions ===\n";
    std::cout << "PDB ID: " << pdbID << "\n";

    // Extract info for all chains
    auto chainInfos = extractChainInfo(pdbID, moleculeData);

    if (chainInfos.empty()) {
        std::cout << "\nNo valid chains found with UniProt mappings.\n";
        std::cout << "Cannot suggest mutations.\n\n";
        return allSuggestions;
    }

    // For each chain, get all pathogenic mutations
    for (const auto& chainInfo : chainInfos) {
        auto chainMutations = getAllMutationsForUniProt(
            chainInfo.uniprotID,
            chainInfo.chainID,
            chainInfo.sequence,
            chainInfo.positions,
            0.7  // Min pathogenicity threshold
        );

        // Add to all suggestions
        allSuggestions.insert(
            allSuggestions.end(),
            chainMutations.begin(),
            chainMutations.end()
        );
    }

    if (allSuggestions.empty()) {
        std::cout << "\nNo pathogenic mutations found in AlphaMissense database.\n";
        std::cout << "This protein may not have well-characterized disease mutations.\n\n";
        return allSuggestions;
    }

    // Sort by pathogenicity score (highest first)
    std::sort(allSuggestions.begin(), allSuggestions.end(),
        [](const MutationSuggestion& a, const MutationSuggestion& b) {
            return a.pathogenicityScore > b.pathogenicityScore;
        }
    );

    // Limit to max suggestions
    if (allSuggestions.size() > maxSuggestions) {
        allSuggestions.resize(maxSuggestions);
    }

    // Display suggestions
    std::cout << "\n=== Top " << allSuggestions.size() << " Suggested Mutations ===\n";
    for (size_t i = 0; i < allSuggestions.size(); ++i) {
        std::cout << (i + 1) << ". " << allSuggestions[i].toString() << "\n";
    }

    std::cout << "\nUse command: mutate <chain> <position> <amino_acid>\n";
    std::cout << "Example: mutate " << allSuggestions[0].chain << " "
        << allSuggestions[0].position << " "
        << allSuggestions[0].suggestedAA << "\n\n";

    return allSuggestions;
}