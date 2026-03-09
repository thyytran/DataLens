#include "UniProtService.h"
#include <iostream>
#include <sstream>
#include <windows.h>
#undef min
#undef max

// Initialize amino acid properties lookup table
const std::map<char, AminoAcidProperties> UniProtService::AA_PROPERTIES = {
	{'A', {"Alanine", "ALA", "A", "neutral", 88.6, 1.8, "nonpolar", 6.00}},
	{'R', {"Arginine", "ARG", "R", "positive", 173.4, -4.5, "polar", 10.76}},
	{'N', {"Asparagine", "ASN", "N", "neutral", 114.1, -3.5, "polar", 5.41}},
	{'D', {"Aspartic acid", "ASP", "D", "negative", 111.1, -3.5, "charged", 2.77}},
	{'C', {"Cysteine", "CYS", "C", "neutral", 108.5, 2.5, "polar", 5.07}},
	{'Q', {"Glutamine", "GLN", "Q", "neutral", 143.8, -3.5, "polar", 5.65}},
	{'E', {"Glutamic acid", "GLU", "E", "negative", 138.4, -3.5, "charged", 3.22}},
	{'G', {"Glycine", "GLY", "G", "neutral", 60.1, -0.4, "nonpolar", 5.97}},
	{'H', {"Histidine", "HIS", "H", "positive (at pH 6)", 153.2, -3.2, "polar", 7.59}},
	{'I', {"Isoleucine", "ILE", "I", "neutral", 166.7, 4.5, "nonpolar", 6.02}},
	{'L', {"Leucine", "LEU", "L", "neutral", 166.7, 3.8, "nonpolar", 5.98}},
	{'K', {"Lysine", "LYS", "K", "positive", 168.6, -3.9, "charged", 9.74}},
	{'M', {"Methionine", "MET", "M", "neutral", 162.9, 1.9, "nonpolar", 5.74}},
	{'F', {"Phenylalanine", "PHE", "F", "neutral", 189.9, 2.8, "nonpolar", 5.48}},
	{'P', {"Proline", "PRO", "P", "neutral", 112.7, -1.6, "nonpolar", 6.30}},
	{'S', {"Serine", "SER", "S", "neutral", 89.0, -0.8, "polar", 5.68}},
	{'T', {"Threonine", "THR", "T", "neutral", 116.1, -0.7, "polar", 5.60}},
	{'W', {"Tryptophan", "TRP", "W", "neutral", 227.8, -0.9, "nonpolar", 5.89}},
	{'Y', {"Tyrosine", "TYR", "Y", "neutral", 193.6, -1.3, "polar", 5.66}},
	{'V', {"Valine", "VAL", "V", "neutral", 140.0, 4.2, "nonpolar", 5.96}}
};

/*
 * Constructor
 */
UniProtService::UniProtService(const std::string& serviceURL)
	: serviceURL(serviceURL), dataLoaded(false) {
}

/*
 * Make HTTP GET request to service
 */
std::string UniProtService::makeRequest(const std::string& endpoint) {
	std::string url = serviceURL + endpoint;

	if (!httpConnection.get(url)) {
		return "";
	}

	return httpConnection.getResponse();
}

/*
 * Check if service is running
 */
bool UniProtService::checkServiceHealth() {
	std::string response = makeRequest("/");

	if (response.empty()) {
		std::cerr << "Error: UniProt service not responding at " << serviceURL << "\n";
		std::cerr << "Make sure the service is running: python3 uniprot_service.py\n";
		return false;
	}

	try {
		json healthData = json::parse(response);
		if (healthData["status"] == "running") {
			std::cout << "UniProt service is running ✓\n";
			return true;
		}
	}
	catch (...) {
		std::cerr << "Error: Invalid response from service\n";
		return false;
	}

	return false;
}

/*
 * Fetch full protein data from service
 */
bool UniProtService::fetchProteinData(const std::string& uniprotID) {
	this->currentUniProtID = uniprotID;

	std::cout << "Fetching data for " << uniprotID << " from UniProt service...\n";

	std::string endpoint = "/protein/" + uniprotID;
	std::string response = makeRequest(endpoint);

	if (response.empty()) {
		std::cerr << "Error: No response from service\n";
		return false;
	}

	try {
		cachedProteinData = json::parse(response);
		dataLoaded = true;

		std::cout << "✓ Loaded: " << cachedProteinData["protein_name"] << "\n";
		std::cout << "  Sequence length: " << cachedProteinData["sequence_length"] << " residues\n";

		return true;
	}
	catch (json::exception& e) {
		std::cerr << "Error parsing response: " << e.what() << "\n";
		return false;
	}
}

/*
 * Get residue-specific information from service
 */
std::optional<ResidueInfo> UniProtService::getResidueInfo(int position) {
	if (currentUniProtID.empty()) {
		std::cerr << "Error: No UniProt ID set. Call fetchProteinData() first.\n";
		return std::nullopt;
	}

	std::string endpoint = "/protein/" + currentUniProtID + "/residue/" + std::to_string(position);
	std::string response = makeRequest(endpoint);

	if (response.empty()) {
		std::cerr << "Error: No response from service for residue " << position << "\n";
		return std::nullopt;
	}

	try {
		json residueData = json::parse(response);

		ResidueInfo info;
		info.position = residueData["position"];
		info.aminoAcid = residueData["amino_acid"].get<std::string>()[0];

		// Parse secondary structure
		info.secondaryStructure.type = residueData["secondary_structure"]["type"];
		std::string range = residueData["secondary_structure"]["range"];
		if (!range.empty()) {
			size_t dashPos = range.find('-');
			if (dashPos != std::string::npos) {
				info.secondaryStructure.start = std::stoi(range.substr(0, dashPos));
				info.secondaryStructure.end = std::stoi(range.substr(dashPos + 1));
			}
		}

		info.isDNABinding = residueData["dna_binding"];

		// Parse variants
		for (const auto& var : residueData["variants"]) {
			Variant v;
			v.position = var["position"];
			v.fromAA = var["from_aa"].get<std::string>()[0];
			v.toAA = var["to_aa"].get<std::string>()[0];
			v.description = var["description"];
			v.disease = var["disease"];
			info.variants.push_back(v);
		}

		// Parse PTMs
		for (const auto& ptm : residueData["ptms"]) {
			PTM p;
			p.position = ptm["position"];
			p.modification = ptm["modification"];
			info.ptms.push_back(p);
		}

		// Get amino acid properties
		info.properties = getAAProperties(info.aminoAcid);

		return info;
	}
	catch (json::exception& e) {
		std::cerr << "Error parsing residue info: " << e.what() << "\n";
		return std::nullopt;
	}
}

/*
 * Get protein name
 */
std::string UniProtService::getProteinName() const {
	if (!dataLoaded) return "";

	try {
		std::string name = cachedProteinData["protein_name"];
		return name != "null" ? name : "";
	}
	catch (...) {
		return "";
	}
}

/*
 * Get protein function
 */
std::string UniProtService::getFunction() const {
	if (!dataLoaded) return "";

	try {
		if (cachedProteinData["function"].is_null()) {
			return "";
		}
		return cachedProteinData["function"];
	}
	catch (...) {
		return "";
	}
}

/*
 * Get protein sequence
 */
std::string UniProtService::getSequence() const {
	if (!dataLoaded) return "";

	try {
		if (cachedProteinData["sequence"].is_null()) {
			return "";
		}
		return cachedProteinData["sequence"];
	}
	catch (...) {
		return "";
	}
}

/*
 * Get sequence length
 */
int UniProtService::getSequenceLength() const {
	if (!dataLoaded) return 0;

	try {
		return cachedProteinData["sequence_length"];
	}
	catch (...) {
		return 0;
	}
}

/*
 * Get current UniProt ID
 */
std::string UniProtService::getCurrentUniProtID() const {
	return currentUniProtID;
}

/*
 * Parse secondary structure from cached data
 */
std::vector<SecondaryStructure> UniProtService::parseSecondaryStructure(const json& data) {
	std::vector<SecondaryStructure> result;

	try {
		for (const auto& ss : data["secondary_structure"]) {
			SecondaryStructure s;
			s.type = ss["type"];
			s.start = ss["start"];
			s.end = ss["end"];
			result.push_back(s);
		}
	}
	catch (...) {
		// Return empty if parsing fails
	}

	return result;
}

/*
 * Get secondary structure type for a position
 */
std::string UniProtService::getSecondaryStructureType(int position) {
	if (!dataLoaded) return "Unknown";

	auto ssList = parseSecondaryStructure(cachedProteinData);

	for (const auto& ss : ssList) {
		if (ss.start <= position && position <= ss.end) {
			return ss.type;
		}
	}

	return "Loop/Coil";
}

/*
 * Check if position is in DNA binding region
 */
bool UniProtService::isDNABinding(int position) {
	if (!dataLoaded) return false;

	try {
		for (const auto& binding : cachedProteinData["dna_binding"]) {
			int start = binding["start"];
			int end = binding["end"];
			if (start <= position && position <= end) {
				return true;
			}
		}
	}
	catch (...) {
		// Return false if parsing fails
	}

	return false;
}

/*
 * Get amino acid properties by one-letter code
 */
AminoAcidProperties UniProtService::getAAProperties(char oneLetterCode) {
	auto it = AA_PROPERTIES.find(toupper(oneLetterCode));
	if (it != AA_PROPERTIES.end()) {
		return it->second;
	}

	// Return default/unknown if not found
	return { "Unknown", "???", "X", "unknown", 0.0, 0.0, "unknown", 7.0 };
}

/*
 * Print residue information (for debugging/display)
 */
void UniProtService::printResidueInfo(const ResidueInfo& info) const {
	std::cout << "\n=== RESIDUE INFORMATION ===\n";
	std::cout << "Position: " << info.position << "\n";
	std::cout << "Residue: " << info.properties.fullName
		<< " (" << info.properties.threeLetterCode << ", "
		<< info.aminoAcid << ")\n";

	std::cout << "\n--- Biochemical Properties ---\n";
	std::cout << "Charge: " << info.properties.charge << "\n";
	std::cout << "Size: " << info.properties.size << " Ų\n";
	std::cout << "Hydrophobicity: " << info.properties.hydrophobicity
		<< " (" << info.properties.category << ")\n";
	std::cout << "Isoelectric point: " << info.properties.pI << "\n";

	std::cout << "\n--- Structural Context ---\n";
	std::cout << "Secondary structure: " << info.secondaryStructure.type;
	if (!info.secondaryStructure.type.empty() && info.secondaryStructure.type != "Loop/Coil") {
		std::cout << " (residues " << info.secondaryStructure.start
			<< "-" << info.secondaryStructure.end << ")";
	}
	std::cout << "\n";
	std::cout << "DNA binding region: " << (info.isDNABinding ? "Yes" : "No") << "\n";

	if (!info.variants.empty()) {
		std::cout << "\n--- Known Variants at This Position ---\n";
		for (const auto& var : info.variants) {
			std::cout << "  " << var.fromAA << var.position << var.toAA
				<< " → " << var.disease << "\n";
			std::cout << "  " << var.description << "\n";
		}
	}
	else {
		std::cout << "\n--- Known Variants ---\n";
		std::cout << "No known disease variants at this position\n";
	}

	if (!info.ptms.empty()) {
		std::cout << "\n--- Post-Translational Modifications ---\n";
		for (size_t i = 0; i < std::min(info.ptms.size(), size_t(3)); ++i) {
			std::cout << "  • " << info.ptms[i].modification << "\n";
		}
		if (info.ptms.size() > 3) {
			std::cout << "  ... and " << (info.ptms.size() - 3) << " more\n";
		}
	}

	std::cout << "=========================\n\n";
}