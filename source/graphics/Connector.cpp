#include "graphics/Connector.h"

/*
 *  Constructor for the Connector class.
 *  This constructor creates a connector (cylinder) between two atoms.
 *  @param atom1             A pointer to the first atom.
 *  @param atom2             A pointer to the second atom.
 *  @param radius            The radius of the connector (cylinder).
 *  @param topColor          A pointer to the color of the top of the connector.
 *  @param bottomColor       A pointer to the color of the bottom of the connector.
 *  @param point1            The 3D coordinates of the first atom.
 *  @param point2            The 3D coordinates of the second atom.
 */
Connector::Connector(
	const Atom *atom1, 
	const Atom *atom2,
	float radius,
	const Color *topColor, 
	const Color *bottomColor, 
	const Vec3 &point1, const Vec3 &point2) :

	atom1(atom1), // Initialize the atom1 member variable
	atom2(atom2), // Initialize the atom2 member variable
	radius(radius), // Initialize the radius member variable
	topColor(*topColor), // Initialize the topColor member variable by copying the value pointed to by topColor
	bottomColor(*bottomColor) // Initialize the bottomColor member variable by copying the value pointed to by bottomColor
{
	// Calculate the length of the connector (distance between the two atoms)
	length = std::sqrt(
		(point1.getX() - point2.getX()) * (point1.getX() - point2.getX()) +
		(point1.getY() - point2.getY()) * (point1.getY() - point2.getY()) +
		(point1.getZ() - point2.getZ()) * (point1.getZ() - point2.getZ())
	);

	Vec3 midPoint(
		(point1.getX() + point2.getX()) / 2,
		(point1.getY() + point2.getY()) / 2,
		(point1.getZ() + point2.getZ()) / 2
	);

	// Calculate the rotation axis and angle to align the connector with the two atoms
	Vec3 defaultDirection(0.0f, 1.0f, 0.0f); // Define a default direction vector (pointing up along the y-axis)

	// Calculate the vector pointing from atom2 to atom1
	Vec3 pointDifference = point1 - point2;

	// Calculate the rotation axis by taking the cross product of the default direction and the point difference
	Vec3 rotationAxis = defaultDirection.cross(pointDifference);

	float rotationAngle = std::acos(  // Calculate the rotation angle using the arccosine function
		//-1 to +1: when +1 no rotation is needed
		defaultDirection.dot(pointDifference) / pointDifference.mag()
	);

	scaleMatrix = MathUtils::MatGen::scale<float, 4>(Vec3(radius, length, radius));
	rotationMatrix = MathUtils::MatGen::rotationAboutAxis<float>(rotationAxis, rotationAngle);
	translationMatrix = MathUtils::MatGen::translation<float, 4>(midPoint);
}

/*
 *  Renders the connector (cylinder).
 *  @param shader              A pointer to the shader program to use for rendering.
 *  @param modelRotationMatrix A matrix representing the rotation of the entire model.
 *  @param connectorTemplate   A pointer to the ConnectorTemplate object containing the vertex data for the cylinder.
 */
void Connector::render(
	const Shader *shader,
	const Mat4 &modelRotationMatrix,
	const ConnectorTemplate *connectorTemplate
) {
	Mat4 modelMatrix = scaleMatrix * rotationMatrix * translationMatrix * modelRotationMatrix;

	shader->setModelMatrix(modelMatrix);
	shader->setNormalMatrix(MathUtils::MatGen::normal<float>(modelMatrix));
	
	shader->setVec3("topColor", topColor.r, topColor.g, topColor.b);
	shader->setVec3("bottomColor", bottomColor.r, bottomColor.g, bottomColor.b);

	glDrawArrays(
		GL_TRIANGLES,
		0,
		GLsizei(
			connectorTemplate->getVerticesLength() / ConnectorTemplate::VERTICES_PER_POINT
		)
	);
}

/*
 *  Sets the radius of the connector (cylinder).
 *  @param radius The new radius of the connector.
 */
void Connector::setRadius(float radius) {
	this->radius = radius;
	scaleMatrix = MathUtils::MatGen::scale<float, 4>(Vec3(radius, length, radius));
}

/*
 *  Sets the colors of the top and bottom of the connector (cylinder).
 *  @param topColor    A pointer to the new color of the top of the connector. If nullptr, the top color is not changed.
 *  @param bottomColor A pointer to the new color of the bottom of the connector. If nullptr, the bottom color is not changed.
 */
void Connector::setColors(const Color *topColor, const Color *bottomColor) {
	if (topColor) {
		this->topColor = *topColor;
	}
	if (bottomColor) {
		this->bottomColor = *bottomColor;
	}
}
