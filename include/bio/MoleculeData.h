#pragma once

#include <vector>
#include <iostream>

#include "Residue.h"
#include "Chain.h"
#include "Helix.h"
#include "Sheet.h"
#include "DisulfideBond.h"
#include "Atom.h"

// Class representing the data for a molecule (e.g., a protein)
class MoleculeData {
public:

	// Vector to store the sequence of [Residue, Chain, Helix, Sheet, DisulfideBond, Atom] in the molecule
	std::vector<Residue> sequence;
	std::vector<Chain> chains;
	std::vector<Helix> helices;
	std::vector<Sheet> sheets;
	std::vector<DisulfideBond> disulfideBonds;
	std::vector<Atom> atoms;

	//virtual const std::vector<Atom*>& getAtoms() const = 0;  // Pure virtual

	void printSequence();
	void printHelices();
	void printSheets();
	void printDisulfideBonds();
	void printAtoms();
	void printChains();
};
