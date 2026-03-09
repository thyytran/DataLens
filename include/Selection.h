#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <optional>

#include "Parser.h"
#include "bio/Atom.h"

class Selection {
public:
	std::optional<int> residue;
	std::optional<std::pair<std::optional<int>, std::optional<int>>> residueRange;	
	std::optional<std::string> element;
	std::optional<char> chain;

	Selection();

	void reset();
	void print() const;
	void parseQuery(const std::vector<std::string> &query);

	bool isMatch(const Atom *atom, bool reversed = false) const;
};
