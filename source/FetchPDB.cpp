#include "FetchPDB.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <iomanip>

using json = nlohmann::json;

const std::string FetchPDB::API_BASE_URL = "http://localhost:8000";

std::string FetchPDB::buildPDBUrl(const std::string& pdb_id) {
    if (pdb_id.size() == 4) {
        return "http://files.rcsb.org/view/" + pdb_id + ".pdb";
    }
    return pdb_id;
}

std::string FetchPDB::buildAPIUrl(const std::string& pdb_id) {
    return API_BASE_URL + "/api/structure/" + pdb_id;
}

PDBStructureData FetchPDB::parseStructureResponse(const std::string& json_response) {
    PDBStructureData data;

    try {
        json response = json::parse(json_response);

        data.pdb_id = response["pdb_id"];
        data.total_mutations = response["total_mutations"];
        data.has_data = true;

        // Parse chains
        if (response.contains("chains")) {
            for (const auto& chain_json : response["chains"]) {
                PDBChainMapping chain;
                chain.chain_id = chain_json["chain_id"];
                chain.uniprot_id = chain_json["uniprot_id"];
                chain.pdb_start = chain_json["pdb_start"];
                chain.pdb_end = chain_json["pdb_end"];
                chain.uniprot_start = chain_json["uniprot_start"];
                chain.uniprot_end = chain_json["uniprot_end"];
                chain.coverage_length = chain_json["coverage_length"];
                data.chains.push_back(chain);
            }
        }

        // Parse summaries
        if (response.contains("chain_summaries")) {
            for (const auto& summary_json : response["chain_summaries"]) {
                PDBMutationSummary summary;
                summary.chain_id = summary_json["chain_id"];
                summary.uniprot_id = summary_json["uniprot_id"];
                summary.uniprot_start = summary_json["uniprot_start"];
                summary.uniprot_end = summary_json["uniprot_end"];
                summary.positions_with_data = summary_json["positions_with_data"];
                summary.total_mutations = summary_json["total_mutations"];
                summary.avg_pathogenicity = summary_json["avg_pathogenicity"];
                summary.pathogenic_count = summary_json["pathogenic_count"];
                summary.benign_count = summary_json["benign_count"];
                summary.ambiguous_count = summary_json["ambiguous_count"];
                data.summaries.push_back(summary);
            }
        }

    }
    catch (const json::exception& e) {
        std::cerr << "JSON parsing error: " << e.what() << std::endl;
        data.has_data = false;
    }

    return data;
}

MoleculeData* FetchPDB::fetchPDBStructure(const std::string& pdb_id) {
    std::string url = buildPDBUrl(pdb_id);
    std::cout << "Downloading PDB structure from: " << url << std::endl;

    try {
        MoleculeData* moleculeData = new PDBFile(url);
        return moleculeData;
    }
    catch (const std::exception& e) {
        std::cerr << "Error loading PDB: " << e.what() << std::endl;
        return nullptr;
    }
}

PDBStructureData FetchPDB::fetchStructureData(const std::string& pdb_id) {
    HTTPConnection conn;
    std::string url = buildAPIUrl(pdb_id);

    std::cout << "Fetching mutation data from API..." << std::endl;

    if (!conn.get(url)) {
        std::cerr << "Failed to fetch structure data from API" << std::endl;
        return PDBStructureData();
    }

    return parseStructureResponse(conn.response);
}

void FetchPDB::displayStructureData(const PDBStructureData& data) {
    if (!data.has_data) {
        std::cout << "\nNo mutation data available for this structure." << std::endl;
        return;
    }

    std::cout << "\n";
    std::cout << "==========================================" << std::endl;
    std::cout << "  PDB STRUCTURE: " << data.pdb_id << std::endl;
    std::cout << "==========================================" << std::endl;

    // Display chains
    std::cout << "\nCHAINS (" << data.chains.size() << "):" << std::endl;
    std::cout << "==========================================" << std::endl;

    for (const auto& chain : data.chains) {
        std::cout << "Chain " << chain.chain_id << " | "
            << std::setw(10) << std::left << chain.uniprot_id << " | "
            << "PDB: " << std::setw(4) << chain.pdb_start << "-"
            << std::setw(4) << std::left << chain.pdb_end << " | "
            << "UniProt: " << std::setw(4) << chain.uniprot_start << "-"
            << std::setw(4) << chain.uniprot_end
            << " (" << chain.coverage_length << " residues)" << std::endl;
    }

    // Display 


    std::cout << "\nMUTATION ANALYSIS:" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Total mutations available: " << data.total_mutations << std::endl;
    std::cout << std::endl;

    for (const auto& summary : data.summaries) {
        std::cout << "Chain " << summary.chain_id << " (" << summary.uniprot_id << "):" << std::endl;
        std::cout << "  Positions covered: " << summary.positions_with_data << std::endl;
        std::cout << "  Total mutations: " << summary.total_mutations << std::endl;
        std::cout << "  Avg pathogenicity: " << std::fixed << std::setprecision(3)
            << summary.avg_pathogenicity << std::endl;
        std::cout << "  Pathogenic: " << summary.pathogenic_count
            << " | Benign: " << summary.benign_count
            << " | Ambiguous: " << summary.ambiguous_count << std::endl;
        std::cout << std::endl;
    }

    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
}

bool FetchPDB::fetchComplete(
    const std::string& pdb_id,
    MoleculeData*& moleculeData,
    PDBStructureData& structureData
) {
    std::cout << "\n[1/3] Fetching PDB structure..." << std::endl;

    moleculeData = fetchPDBStructure(pdb_id);
    if (moleculeData == nullptr) {
        std::cerr << "Failed to load PDB structure" << std::endl;
        return false;
    }

    std::cout << "[2/3] Fetching mutation data..." << std::endl;

    structureData = fetchStructureData(pdb_id);

    std::cout << "[3/3] Loading complete!" << std::endl;

    displayStructureData(structureData);

    if (structureData.has_data && structureData.total_mutations > 0) {
        std::cout << "Structure loaded with "
            << structureData.total_mutations
            << " mutations available for analysis" << std::endl;
        return true;
    }
    else {
        std::cout << "Structure loaded but no mutation data available" << std::endl;
        return true;
    }
}