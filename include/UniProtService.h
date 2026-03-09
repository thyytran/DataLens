#ifndef UNIPROT_SERVICE_H
#define UNIPROT_SERVICE_H

#include <string>
#include <vector>
#include <map>
#include <optional>
#include "nlohmann/json.hpp"
#include "HTTPConnection.h"

using json = nlohmann::json;

// Amino acid property structure
struct AminoAcidProperties {
	std::string fullName;
	std::string threeLetterCode;
	std::string oneLetterCode;
	std::string charge;       // "positive", "negative", "neutral"
	double size;              // Volume in Ų
	double hydrophobicity;    // Kyte-Doolittle scale
	std::string category;     // "polar", "nonpolar", "charged"
	double pI;                // Isoelectric point
};

// Secondary structure element
struct SecondaryStructure {
	std::string type;  // "Helix", "Beta strand", "Turn", "Loop/Coil"
	int start;
	int end;
};

// Known variant
struct Variant {
	int position;
	char fromAA;
	char toAA;
	std::string description;
	std::string disease;
};

// Post-translational modification
struct PTM {
	int position;
	std::string modification;
};

// Residue information structure
struct ResidueInfo {
	int position;
	char aminoAcid;
	SecondaryStructure secondaryStructure;
	bool isDNABinding;
	std::vector<Variant> variants;
	std::vector<PTM> ptms;
	AminoAcidProperties properties;
};

// Main UniProt service client class
class UniProtService {
private:
	std::string serviceURL;     // e.g., "http://localhost:8000"
	std::string currentUniProtID;
	json cachedProteinData;
	bool dataLoaded;

	// HTTP connection for service requests
	HTTPConnection httpConnection;

	// Static amino acid properties lookup table
	static const std::map<char, AminoAcidProperties> AA_PROPERTIES;

	// Make HTTP GET request to service
	std::string makeRequest(const std::string& endpoint);

	// Parse JSON responses
	std::vector<SecondaryStructure> parseSecondaryStructure(const json& data);
	std::vector<Variant> parseVariants(const json& data);
	std::vector<PTM> parsePTMs(const json& data);

public:
	// Constructor with default service URL
	UniProtService(const std::string& serviceURL = "http://localhost:8000");

	// Fetch full protein data from service
	bool fetchProteinData(const std::string& uniprotID);

	// Get residue-specific information
	std::optional<ResidueInfo> getResidueInfo(int position);

	// Get protein name
	std::string getProteinName() const;

	// Get protein function
	std::string getFunction() const;

	// Get full sequence
	std::string getSequence() const;

	// Get sequence length
	int getSequenceLength() const;

	// Get current UniProt ID
	std::string getCurrentUniProtID() const;

	// Get amino acid properties by one-letter code (static)
	static AminoAcidProperties getAAProperties(char oneLetterCode);

	// Get secondary structure type for a position
	std::string getSecondaryStructureType(int position);

	// Check if position is in DNA binding region
	bool isDNABinding(int position);

	// Print residue info (for debugging/console display)
	void printResidueInfo(const ResidueInfo& info) const;

	// Check if service is running
	bool checkServiceHealth();
};

#endif // UNIPROT_SERVICE_H