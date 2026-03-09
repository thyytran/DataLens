#include "bio/Protein.h"

// Default constructor
Protein::Protein() {}

/*
 *  Constructor that takes a string as input and creates a protein sequence from it.
 *  @param sequence A string representing the amino acid sequence (e.g., "AlaGlySer").
 */
Protein::Protein(std::string sequence) {
	std::stringstream sequenceStream(sequence);
	std::string abbr;
	while (sequenceStream >> abbr) {
		const AminoAcid *aminoAcid = AminoAcid::get(abbr);
		if (aminoAcid) {
			this->sequence.push_back(aminoAcid);
		}
		else {
			std::cout << "INFO > Adding amino acid \"" << abbr << "\" to dictionary. " <<
				"Consider creating a dictionary entry yourself.\n\n";
			const AminoAcid *newAminoAcid = AminoAcid::set(abbr, "", "", "", "");
			if (newAminoAcid) {
				this->sequence.push_back(newAminoAcid);
			}
		}
	}
}

/*
 *  Constructor that takes a vector of AminoAcid pointers as input.
 *  @param sequence A vector of pointers to const AminoAcid objects.
 */
Protein::Protein(const std::vector<const AminoAcid*> &sequence) {
	for (size_t i = 0; i < sequence.size(); ++i) {
		this->sequence.push_back(sequence[i]);
	}
}

/*
 *  Constructor that takes a MoleculeData object as input and extracts the protein sequence from it.
 *  @param moleculeData A pointer to a const MoleculeData object containing the molecular data.
 */
Protein::Protein(const MoleculeData *moleculeData) {
	for (size_t i = 0; i < moleculeData->sequence.size(); ++i) {
		std::string name = moleculeData->sequence[i].name;
		const AminoAcid *aminoAcid = AminoAcid::get(name);
		if (aminoAcid) {
			sequence.push_back(aminoAcid);
		}
		else {
			std::cout << "INFO > Adding amino acid \"" << name << "\" to dictionary. " <<
				"Consider creating a dictionary entry yourself.\n\n";
			const AminoAcid *newAminoAcid = AminoAcid::set(name, "", "", "", "");
			if (newAminoAcid) {
				this->sequence.push_back(newAminoAcid);
			}
		}
	}
}

/*
 *  Clears the sequence vector, effectively resetting the protein.
 */
void Protein::reset() {
	sequence.clear();
}

/*
 *  Overloads the << operator to allow printing a Protein object to an output stream.
 *  @param os The output stream to print to.
 *  @param protein The Protein object to print.
 *  @return A reference to the output stream.
 */
std::ostream &operator<<(std::ostream &os, Protein &protein) {
	if (protein.sequence.size() == 0) {
		os << "Empty";
	}
	else {
		for (size_t i = 0; i < protein.sequence.size(); ++i) {
			os << protein.sequence[i]->abbr3 << " ";
		}
	}
	return os;
}
