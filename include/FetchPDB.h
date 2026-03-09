#ifndef DATALENS_FETCH_PDB_H
#define DATALENS_FETCH_PDB_H

#include <string>
#include <vector>
#include "HTTPConnection.h"
#include "bio/PDBFile.h"
#include "bio/MoleculeData.h"

// Data structures for PDB mapping (renamed to avoid conflicts)
struct PDBChainMapping {
    std::string chain_id;
    std::string uniprot_id;
    std::string pdb_start;
    std::string pdb_end;
    int uniprot_start;
    int uniprot_end;
    int coverage_length;

    PDBChainMapping() : uniprot_start(0), uniprot_end(0), coverage_length(0) {}
};

struct PDBMutationSummary {
    std::string chain_id;
    std::string uniprot_id;
    int uniprot_start;
    int uniprot_end;
    int positions_with_data;
    int total_mutations;
    float avg_pathogenicity;
    int pathogenic_count;
    int benign_count;
    int ambiguous_count;

    PDBMutationSummary() : uniprot_start(0), uniprot_end(0), positions_with_data(0),
        total_mutations(0), avg_pathogenicity(0.0f),
        pathogenic_count(0), benign_count(0), ambiguous_count(0) {
    }
};

struct PDBStructureData {
    std::string pdb_id;
    std::vector<PDBChainMapping> chains;
    std::vector<PDBMutationSummary> summaries;
    int total_mutations;
    bool has_data;

    PDBStructureData() : total_mutations(0), has_data(false) {}
};

class FetchPDB {
private:
    static const std::string API_BASE_URL;

    static std::string buildPDBUrl(const std::string& pdb_id);
    static std::string buildAPIUrl(const std::string& pdb_id);
    static PDBStructureData parseStructureResponse(const std::string& json_response);

public:
    static MoleculeData* fetchPDBStructure(const std::string& pdb_id);
    static PDBStructureData fetchStructureData(const std::string& pdb_id);
    static void displayStructureData(const PDBStructureData& data);
    static bool fetchComplete(
        const std::string& pdb_id,
        MoleculeData*& moleculeData,
        PDBStructureData& structureData
    );
};

#endif // DATALENS_FETCH_PDB_H