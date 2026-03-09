#pragma once

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "HTTPConnection.h"
#include "MoleculeData.h"
#include "AminoAcid.h"
#include "Helix.h"
#include "Sheet.h"
#include "DisulfideBond.h"
#include "Atom.h"
#include "Parser.h"
#include "math/Vec.h"

class PDBFile : public MoleculeData {
private:
    std::string pdb_id;
    void parseFromString(const std::string& pdbText);

public:
	PDBFile(const std::string &url);
    PDBFile(const std::string& source, bool isRawContent);

    ~PDBFile();


    // Optional: add these too
    size_t getAtomCount() const { return atoms.size(); }
    const std::string& getPDBID() const { return pdb_id; }
    // const std::vector<Atom*>& getAtoms() const override { return atoms; }  // ADD THIS

};
