#include "graphics/Color.h"

// Define constant Color objects for different helix types
// const Color Color::HELIX_COLOR_ALPHA(1.0f, 0.0f, 0.5f);
//const Color Color::HELIX_COLOR_3_10(0.63f, 0.0f, 0.5f);
//const Color Color::HELIX_COLOR_PI(0.38f, 0.0f, 0.5f);
//const Color Color::HELIX_COLOR_UNDEFINED(0.7f, 0.7f, 0.7f);
//const Color Color::SHEET_COLOR(1.0f, 0.78f, 0.0f);

const Color Color::HELIX_COLOR_ALPHA(0.95f, 0.3f, 0.6f);      // Vibrant rose pink
const Color Color::HELIX_COLOR_3_10(0.7f, 0.2f, 0.75f);       // Deep violet
const Color Color::HELIX_COLOR_PI(0.5f, 0.15f, 0.7f);         // Rich indigo
const Color Color::HELIX_COLOR_UNDEFINED(0.65f, 0.68f, 0.7f); // Cool silver-blue
const Color Color::SHEET_COLOR(1.0f, 0.75f, 0.15f);           // Warm amber gold
/*
 *  Default constructor for the Color class.
 *  This constructor does not initialize the r, g, and b members.
 */
Color::Color() {}

/*
 *  Constructor for the Color class.
 *  This constructor initializes the r, g, and b members with the provided values.
 *  @param r The red component of the color.
 *  @param g The green component of the color.
 *  @param b The blue component of the color.
 */
Color::Color(float r, float g, float b) :
	r(r), g(g), b(b) {}

/*
 *  Creates a Color object from byte values.
 *  This function takes three unsigned char values (r, g, b) and converts them to
 *  floating-point values between 0.0f and 1.0f by dividing them by 255.0f.
 *  @param r The red component of the color (0-255).
 *  @param g The green component of the color (0-255).
 *  @param b The blue component of the color (0-255).
 *  @return A Color object with the converted r, g, and b values.
 */
Color Color::fromByte(unsigned char r, unsigned char g, unsigned char b) {
	return Color(r / 255.0f, g / 255.0f, b / 255.0f);
}

/*
 *  Creates a Color object from a color name.
 *  This function takes a string representing a color name and returns a Color object
 *  with the corresponding r, g, and b values.
 *  @param name The name of the color.
 *  @return A Color object with the corresponding r, g, and b values.
 *  @throws std::invalid_argument if the color name is invalid.
 */
Color Color::fromName(const std::string &name) {
	/*
	if (name == "red") {
		return Color(1.0f, 0.0f, 0.0f);
	}
	if (name == "green") {
		return Color(0.0f, 1.0f, 0.0f);
	}
	if (name == "blue") {
		return Color(0.0f, 0.0f, 1.0f);
	}
	if (name == "orange") {
		return Color(1.0f, 0.5f, 0.0f);
	}
	if (name == "yellow") {
		return Color(1.0f, 1.0f, 0.0f);
	}
	if (name == "purple") {
		return Color(1.0f, 0.0f, 1.0f);
	}
	if (name == "white") {
		return Color(1.0f, 1.0f, 1.0f);
	}
	if (name == "light-gray") {
		return Color(0.7f, 0.7f, 0.7f);
	}
	if (name == "dark-gray") {
		return Color(0.3f, 0.3f, 0.3f);
	}
	if (name == "brown") {
		return Color(0.57f, 0.36f, 0.19f);
	}
	if (name == "black") {
		return Color(0.0f, 0.0f, 0.0f);
	}
	*/

	if (name == "red") {
		// Vibrant cherry red
		return Color(0.9f, 0.2f, 0.25f);
	}
	if (name == "green") {
		// Rich emerald green
		return Color(0.2f, 0.8f, 0.3f);
	}
	if (name == "blue") {
		// Deep azure blue
		return Color(0.15f, 0.45f, 0.9f);
	}
	if (name == "orange") {
		// Warm sunset orange
		return Color(1.0f, 0.55f, 0.2f);
	}
	if (name == "yellow") {
		// Bright sunflower yellow
		return Color(0.95f, 0.9f, 0.2f);
	}
	if (name == "purple") {
		// Royal amethyst purple
		return Color(0.7f, 0.25f, 0.85f);
	}
	if (name == "white") {
		// Soft cream white
		return Color(0.98f, 0.98f, 0.96f);
	}
	if (name == "light-gray") {
		// Cool light silver
		return Color(0.75f, 0.77f, 0.78f);
	}
	if (name == "dark-gray") {
		// Warm slate gray
		return Color(0.35f, 0.36f, 0.38f);
	}
	if (name == "brown") {
		// Rich mahogany brown
		return Color(0.6f, 0.35f, 0.25f);
	}
	if (name == "black") {
		// Deep charcoal black
		return Color(0.08f, 0.08f, 0.1f);
	}

	throw std::invalid_argument("Invalid color name");
}

/*
 *  Creates a Color object from an element symbol.
 *  This function takes a string representing an element symbol and returns a Color object
 *  with the corresponding r, g, and b values based on the element.
 *  @param element The element symbol.
 *  @return A Color object with the corresponding r, g, and b values.
 *  @throws std::invalid_argument if the element symbol is unknown.
 */
Color Color::fromElement(const std::string &element) {
	/*
	if (element == "C") {
		return Color(0.39f, 0.39f, 0.39f);
	}
	if (element == "H") {
		return Color(1.0f, 1.0f, 1.0f);
	}
	if (element == "N") {
		return Color(0.0f, 0.0f, 1.0f);
	}
	if (element == "O") {
		return Color(1.0f, 0.0f, 0.0f);
	}
	if (element == "P") {
		return Color(1.0f, 0.5f, 0.0f);
	}
	if (element == "S") {
		return Color(1.0f, 1.0f, 0.0f);
	}
	*/

	if (element == "C") {
		// Sophisticated charcoal with slight warm tint
		return Color(0.25f, 0.24f, 0.22f);
	}
	if (element == "H") {
		// Soft pearl white with hint of blue
		return Color(0.95f, 0.97f, 1.0f);
	}
	if (element == "N") {
		// Deep electric blue with vibrancy
		return Color(0.2f, 0.4f, 0.95f);
	}
	if (element == "O") {
		// Rich crimson red
		return Color(0.85f, 0.15f, 0.2f);
	}
	if (element == "P") {
		// Vibrant tangerine orange
		return Color(1.0f, 0.6f, 0.15f);
	}
	if (element == "S") {
		// Golden yellow with depth
		return Color(0.95f, 0.85f, 0.25f);
	}

	throw std::invalid_argument("Err > Unknown color for element: " + element);
}

/*
 *  Creates a Color object based on the structure of the atom.
 *  This function determines the color of an atom based on its secondary structure
 *  (helix or sheet) within the molecule data.
 *  @param atom A pointer to the Atom object.
 *  @param moleculeData A pointer to the MoleculeData object.
 *  @return A Color object based on the atom's structure.
 */
Color Color::fromStructure(const Atom *atom, const MoleculeData *moleculeData) {
	//Check if atom is in a helix
	for (size_t i = 0; i < moleculeData->helices.size(); ++i) {
		if (atom->chain == moleculeData->helices[i].chain &&
			atom->residueNum >= moleculeData->helices[i].residueStart &&
			atom->residueNum <= moleculeData->helices[i].residueEnd) {

			switch (moleculeData->helices[i].type) {
			//Alpha
			case 1:
			case 6:
				return HELIX_COLOR_ALPHA;
			//3/10
			case 5:
				return HELIX_COLOR_3_10;
			//Pi
			case 3:
				return HELIX_COLOR_PI;
			default:
				return HELIX_COLOR_UNDEFINED;
			}

			break;
		}
	}

	//Check if atom is in a sheet
	for (size_t i = 0; i < moleculeData->sheets.size(); ++i) {
		if (atom->chain == moleculeData->sheets[i].chain &&
			atom->residueNum >= moleculeData->sheets[i].residueStart &&
			atom->residueNum <= moleculeData->sheets[i].residueEnd) {

			return SHEET_COLOR;
		}
	}

	return Color(1.0f, 1.0f, 1.0f);
}

/*
 *  Equality operator for the Color class.
 *  This function compares two Color objects and returns true if their r, g, and b
 *  members are equal, and false otherwise.
 *  @param color The Color object to compare with.
 *  @return True if the two Color objects are equal, false otherwise.
 */
bool Color::operator==(const Color &color) const {
	return r == color.r && g == color.g && b == color.b;
}
