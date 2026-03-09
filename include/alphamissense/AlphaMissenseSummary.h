#pragma once
#include <string>
#include <vector>

struct ResidueAmSummary {
	int pdbResidueNum;
	int uniprotPos;
	float worstScore = -1.0f;
	std::string worstClass;
	std::string worstVariant;
};

struct ChainAmSummary {
	std::string uniprotId;
	std::vector<ResidueAmSummary> residues; // sorted by pdbResidueNum
	bool loaded = false;
};

struct ChainResidueInfo {
	std::string chainId;
	std::string uniprotId;
	int pdbStart;
	int pdbEnd;
	bool hasAmData;
};