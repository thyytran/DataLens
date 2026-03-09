// PDBUniProtMapper.cpp
#include "mapper/PDBToUniProtMapper.h"
#include <sstream>
#include <algorithm>

bool PDBUniProtMapper::loadFromCSV(const std::string& filepath) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Err > Failed to open PDB-UniProt mapping file: "
            << filepath << "\n\n";
        return false;
    }

    std::cout << "Loading PDB-UniProt mappings...\n";

    std::string line;

    // Skip header
    if (!std::getline(file, line)) {
        std::cerr << "Err > Empty mapping file\n\n";
        return false;
    }

    int lineCount = 0;
    int errorCount = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::vector<std::string> fields = Parser::split(line, ',');

        if (fields.size() < 9) {
            errorCount++;
            continue;
        }

        try {
            std::string pdbID = Parser::lowercase(fields[0]);
            char chain = fields[1][0];
            std::string uniprotID = fields[2];

            int pdbBeg = std::stoi(fields[5]);
            int pdbEnd = std::stoi(fields[6]);
            int spBeg = std::stoi(fields[7]);
            int spEnd = std::stoi(fields[8]);

            AlignmentSegment segment;
            segment.pdbStart = pdbBeg;
            segment.pdbEnd = pdbEnd;
            segment.uniprotStart = spBeg;
            segment.uniprotEnd = spEnd;
            segment.uniprotID = uniprotID;

            std::string key = makeKey(pdbID, chain);

            auto it = mappings.find(key);
            if (it == mappings.end()) {
                // Create new mapping
                ChainMapping newMapping;
                newMapping.pdbID = pdbID;
                newMapping.chain = chain;
                newMapping.uniprotID = uniprotID;
                newMapping.segments.push_back(segment);

                mappings[key] = std::move(newMapping);  // Use move
            }
            else {
                // Add to existing
                it->second.segments.push_back(segment);
            }

            lineCount++;

            if (lineCount % 100000 == 0) {
                std::cout << "  Loaded " << lineCount / 1000 << "K entries...\n";
            }
        }
        catch (const std::exception&) {
            errorCount++;
        }
    }

    file.close();
    isLoaded = true;

    std::cout << "Loaded " << lineCount << " entries\n";
    std::cout << "Unique chains: " << mappings.size() << "\n\n";

    return true;
}

std::optional<int> PDBUniProtMapper::pdbToUniProtPosition(
    const std::string& pdbID,
    char chain,
    int pdbPosition
) const {
    if (!isLoaded) return std::nullopt;

    std::string key = makeKey(Parser::lowercase(pdbID), chain);
    auto it = mappings.find(key);

    if (it == mappings.end()) return std::nullopt;

    for (const auto& segment : it->second.segments) {
        if (segment.contains(pdbPosition)) {
            return segment.toUniProtPosition(pdbPosition);
        }
    }

    return std::nullopt;
}

std::optional<std::string> PDBUniProtMapper::getUniProtID(
    const std::string& pdbID,
    char chain
) const {
    if (!isLoaded) return std::nullopt;

    std::string key = makeKey(Parser::lowercase(pdbID), chain);
    auto it = mappings.find(key);

    if (it != mappings.end()) {
        return it->second.uniprotID;
    }

    return std::nullopt;
}

std::optional<ChainMapping> PDBUniProtMapper::getChainMapping(
    const std::string& pdbID,
    char chain
) const {
    if (!isLoaded) return std::nullopt;

    std::string key = makeKey(Parser::lowercase(pdbID), chain);
    auto it = mappings.find(key);

    if (it != mappings.end()) {
        return it->second;
    }

    return std::nullopt;
}

std::vector<ChainMapping> PDBUniProtMapper::getAllChainsForPDB(const std::string& pdbID) const {
    std::vector<ChainMapping> chains;

    if (!isLoaded) return chains;

    std::string lowerPdbID = Parser::lowercase(pdbID);

    // Iterate through all mappings and find matches
    for (const auto& pair : mappings) {
        if (pair.second.pdbID == lowerPdbID) {
            chains.push_back(pair.second);
        }
    }

    return chains;
}

void PDBUniProtMapper::displayPDBInfo(const std::string& pdbID) const {
    if (!isLoaded) {
        std::cout << "Warning: PDB-UniProt mapper not loaded.\n\n";
        return;
    }

    auto chains = getAllChainsForPDB(pdbID);

    if (chains.empty()) {
        std::cout << "Warning: No UniProt mappings found for PDB " << pdbID << "\n";
        std::cout << "This structure may not have sequence information.\n\n";
        return;
    }

    std::cout << "\n================================================\n";
    std::cout << "  PDB-UniProt Chain Mappings: " << Parser::uppercase(pdbID) << "\n";
    std::cout << "================================================\n\n";

    std::cout << "Found " << chains.size() << " chain(s) with UniProt mappings:\n\n";

    for (const auto& chainMapping : chains) {
        std::cout << "------------------------------------------------\n";
        std::cout << "Chain: " << chainMapping.chain << "\n";
        std::cout << "UniProt ID: " << chainMapping.uniprotID << "\n";
        std::cout << "Alignment segments: " << chainMapping.segments.size() << "\n\n";

        for (size_t i = 0; i < chainMapping.segments.size(); ++i) {
            const auto& seg = chainMapping.segments[i];
            std::cout << "  Segment " << (i + 1) << ":\n";
            std::cout << "    PDB positions:     " << seg.pdbStart << " - " << seg.pdbEnd
                << " (" << (seg.pdbEnd - seg.pdbStart + 1) << " residues)\n";
            std::cout << "    UniProt positions: " << seg.uniprotStart << " - " << seg.uniprotEnd
                << " (" << (seg.uniprotEnd - seg.uniprotStart + 1) << " residues)\n";

            int offset = seg.uniprotStart - seg.pdbStart;
            if (offset == 0) {
                std::cout << "    Offset: None (positions match 1:1)\n";
            }
            else if (offset > 0) {
                std::cout << "    Offset: +" << offset << " (PDB pos + " << offset << " = UniProt pos)\n";
            }
            else {
                std::cout << "    Offset: " << offset << " (PDB pos " << offset << " = UniProt pos)\n";
            }
            std::cout << "\n";
        }
    }

    std::cout << "------------------------------------------------\n";
    std::cout << "Tip: Use 'mutate <chain> <position> <amino_acid>'\n";
    std::cout << "Example: mutate " << chains[0].chain << " 10 A\n\n";
}