#include "graphics/ConnectorTemplate.h"

// Define the height of the connector (cylinder)
const float ConnectorTemplate::HEIGHT = 1.0f;

// Define the default radius of the connector (cylinder)
const float ConnectorTemplate::DEFAULT_RADIUS = 0.01f;


/*
 *  Adds a vertex point to the vertex array.
 *  @param point The 3D position of the vertex.
 *  @param normal The normal vector of the vertex.
 *  @param isTop A boolean indicating whether the vertex is on the top or bottom of the cylinder (used for texturing or other effects).
 */
void ConnectorTemplate::addVertexPoint(const Vec3 &point, const Vec3 &normal, bool isTop) {
	vertices[verticesIndex] = point.getX();
	vertices[verticesIndex + 1] = point.getY();
	vertices[verticesIndex + 2] = point.getZ();

	vertices[verticesIndex + 3] = normal.getX();
	vertices[verticesIndex + 4] = normal.getY();
	vertices[verticesIndex + 5] = normal.getZ();

	vertices[verticesIndex + 6] = (float)isTop;

	verticesIndex += VERTICES_PER_POINT;
}

/*
 *  Generates the vertices for a cylinder.
 *  The cylinder is created by dividing it into panels around its circumference.
 */
void ConnectorTemplate::genCylinderVertices() {
	// Calculate the angle between each panel
	float panelAngleIntermediate = 2 * MathUtils::PI / NUM_PANELS;
	
	// Define the Y values for the top, middle, and bottom of the cylinder
	const float Y_VALUES[] = { 0.5f, 0.0f, -0.5f };

	// Iterate over the number of panels
	for (unsigned int i = 0; i < NUM_PANELS; ++i) {

		// Calculate the angle of the current panel
		float panelAngle = panelAngleIntermediate * i;

		// Calculate the angle of the next panel
		float nextPanelAngle = panelAngleIntermediate * (i + 1);

		// Calculate the x and z coordinates of the current panel
		float x = std::cos(panelAngle);
		float z = std::sin(panelAngle);

		// Calculate the x and z coordinates of the next panel
		float nextX = std::cos(nextPanelAngle);
		float nextZ = std::sin(nextPanelAngle);

		// Add top and bottom half triangles
		for (unsigned int j = 0; j <= 1; ++j) {

			// Get the y value for the top of the triangle
			const float TOP_Y = Y_VALUES[j];

			// Get the y value for the bottom of the triangle
			const float BOTTOM_Y = Y_VALUES[j + 1];

			//Top left triangle
			addVertexPoint(Vec3(nextX, TOP_Y, nextZ), Vec3(nextX, 0.0f, nextZ), (bool)j); //Top left
			addVertexPoint(Vec3(nextX, BOTTOM_Y, nextZ), Vec3(nextX, 0.0f, nextZ), (bool)j); //Bottom left
			addVertexPoint(Vec3(x, TOP_Y, z), Vec3(x, 0.0f, z), (bool)j); //Top right

			//Bottom right triangle
			addVertexPoint(Vec3(x, TOP_Y, z), Vec3(x, 0.0f, z), (bool)j); //Top right
			addVertexPoint(Vec3(nextX, BOTTOM_Y, nextZ), Vec3(nextX, 0.0f, nextZ), (bool)j); //Bottom left
			addVertexPoint(Vec3(x, BOTTOM_Y, z), Vec3(x, 0.0f, z), (bool)j); //Bottom right
		}
	}
}


/*
 *  Constructor for the ConnectorTemplate class.
 *  This constructor generates the cylinder vertices, creates the vertex array object (VAO),
 *  creates the vertex buffer object (VBO), binds the VAO, binds the VBO,
 *  and sets up the vertex attribute pointers.
 */
ConnectorTemplate::ConnectorTemplate() {
	
	// Generate the cylinder vertices
	genCylinderVertices();

	glGenVertexArrays(1, &vertexArrayID); // Generate the vertex array object (VAO)

	glGenBuffers(1, &vertexBufferID); // Generate the vertex buffer object (VBO)

	glBindVertexArray(vertexArrayID); // Bind the vertex array object (VAO)


	glBindBuffer(GL_ARRAY_BUFFER, vertexBufferID); 	// Bind the vertex buffer object (VBO)

	// Copy the vertex data to the vertex buffer object (VBO)
	glBufferData(GL_ARRAY_BUFFER, verticesIndex * sizeof(float), vertices, GL_STATIC_DRAW);

	//Positions
	glVertexAttribPointer( // Configure the vertex attribute pointer for vertex positions
		0, 3, GL_FLOAT, GL_FALSE, VERTICES_PER_POINT * sizeof(GLfloat),
		(GLvoid*)0
	);
	glEnableVertexAttribArray(0);  // Enable the vertex attribute array for positions

	//Normals
	glVertexAttribPointer( // Configure the vertex attribute pointer for vertex normals
		1, 3, GL_FLOAT, GL_FALSE, VERTICES_PER_POINT * sizeof(GLfloat),  // Attribute index, size, type, normalized, stride, offset
		(GLvoid*)(3 * sizeof(GLfloat)) // Offset of the normal data in the VBO
	); 
	glEnableVertexAttribArray(1); // Enable the vertex attribute array for normals

	//Position indicators
	glVertexAttribPointer( // Configure the vertex attribute pointer for the isTop flag
		2, 1, GL_FLOAT, GL_FALSE, VERTICES_PER_POINT * sizeof(GLfloat), // Attribute index, size, type, normalized, stride, offset
		(GLvoid*)(6 * sizeof(GLfloat))  // Offset of the isTop data in the VBO
	);
	glEnableVertexAttribArray(2); // Enable the vertex attribute array for the isTop flag
}

/*
 *  Gets a pointer to the vertex array.
 *  @return A pointer to the vertex array.
 */
const float *ConnectorTemplate::getVerticesPtr() const {
	return vertices;
}

/*
 *  Gets the length of the vertex array.
 *  @return The length of the vertex array.
 */
size_t ConnectorTemplate::getVerticesLength() const {
	return verticesIndex;
}

/*
 *  Binds the vertex array object (VAO).
 *  This function must be called before rendering the connector.
 */
void ConnectorTemplate::bind() const {
	glBindVertexArray(vertexArrayID);
}
