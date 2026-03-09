#pragma once

#include <bitset>
#include <string>
#include <cctype>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <iomanip>

class AminoAcid { // Class representing an amino acid
private:
	// Static dictionary to store all created AminoAcid objects, accessible by name (full or abbreviated)
	static std::unordered_map<std::string, std::shared_ptr<const AminoAcid>> dict;
	
	/*
	 *  Bitset to store amino acid properties.
	 *  Bit positions:
	 *  0: isPolar - '1' if the amino acid is polar, '0' otherwise
	 *  1: isBasic - '1' if the amino acid is basic, '0' otherwise
	 *  2: isAcidic - '1' if the amino acid is acidic, '0' otherwise
	 */
	const std::bitset<3> properties;

	/*
	 *  Displays a warning message if an attempt is made to overwrite an existing amino acid in the dictionary.
	 *  @param name The name of the amino acid being overwritten.
	 */	
	static void overwriteWarn(const std::string &name);

public:
	// Full name of the amino acid (e.g., "Alanine")
	const std::string name; 

	// One-letter abbreviation of the amino acid (e.g., "A")
	const std::string abbr1;

	// Three-letter abbreviation of the amino acid (e.g., "Ala")
	const std::string abbr3;

	// Chemical formula of the amino acid (e.g., "C3H7NO2")
	const std::string formula;

	/*
	 *  Constructor for the AminoAcid class.
	 *  @param name       The full name of the amino acid.
	 *  @param abbr1      The one-letter abbreviation of the amino acid.
	 *  @param abbr3      The three-letter abbreviation of the amino acid.
	 *  @param formula    The chemical formula of the amino acid.
	 *  @param properties A string of three bits ('0'/'1') representing the amino acid's properties (polar, basic, acidic).
	 */
	AminoAcid(
		std::string name, std::string abbr1, std::string abbr3,
		std::string formula, std::string properties
	);

	bool isPolar() const;
	bool isBasic() const;
	bool isAcidic() const;
	bool containsSulfur() const;

	/*
	 *  Retrieves an AminoAcid object from the dictionary by name (full or abbreviated).
	 *  @param name The name (full or abbreviated) of the amino acid to retrieve.
	 *  @return A pointer to the AminoAcid object if found, nullptr otherwise.
	 */
	static const AminoAcid *get(const std::string &name);

	/*
	 *  Creates a new AminoAcid object, adds it to the dictionary, and returns a pointer to it.
	 *  @param name       The full name of the amino acid.
	 *  @param abbr1      The one-letter abbreviation of the amino acid.
	 *  @param abbr3      The three-letter abbreviation of the amino acid.
	 *  @param formula    The chemical formula of the amino acid.
	 *  @param properties A string of three bits ('0'/'1') representing the amino acid's properties (polar, basic, acidic).
	 *  @return A pointer to the newly created AminoAcid object.
	 */
	static const AminoAcid *set(
		std::string name, std::string abbr1, std::string abbr3,
		std::string formula, std::string properties
	);

	/*
	 *  Removes an amino acid from the dictionary using an AminoAcid pointer.
	 *  @param aminoAcid A pointer to the AminoAcid object to remove.
	 */
	static void erase(const AminoAcid *aminoAcid);

	/*
	 *  Removes an amino acid from the dictionary using its name.
	 *  @param name The name of the amino acid to remove.
	 */
	static void erase(const std::string &name);

	/*
	 *  Displays the contents of the amino acid dictionary to the console.
	 *  This is primarily for debugging purposes.
	 */
	static void showDict();
};
