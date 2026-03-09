#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "HTTPConnection.h"
#include "graphics/SphereTemplate.h"
#include "graphics/Model.h"
#include "AminoAcid.h"
#include "Atom.h"
#include "MoleculeData.h"

class Protein {
public:
	Protein();
	Protein(std::string sequence);
	Protein(const std::vector<const AminoAcid*> &sequence);
	Protein(const MoleculeData *moleculeData);
	void reset();

	std::vector<const AminoAcid*> sequence;
	friend std::ostream &operator<<(std::ostream &os, Protein &protein);
};
