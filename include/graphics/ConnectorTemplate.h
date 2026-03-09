#pragma once

#include <vector>

#include <GL/glew.h>

#include "math/Vec.h"
#include "math/MathUtils.h"

// Class representing a template for creating connector (cylinder) geometry
class ConnectorTemplate {
private:

	// Number of panels (segments) used to create the cylinder around its circumference
	static const unsigned int NUM_PANELS = 14;

	// Height of the cylinder
	static const float HEIGHT;

public:
	// Number of floating-point values per vertex (3 for position, 3 for normal, 1 for side indicator)
	static const unsigned int VERTICES_PER_POINT = 7;


	// Total number of vertices required to create the cylinder
	// Calculated as: 2 (top and bottom halves) * 3 (vertices per triangle) * VERTICES_PER_POINT * NUM_PANELS * 2 (two triangles per panel)
	static const unsigned int NUM_VERTICES = 2 * 3 * VERTICES_PER_POINT * NUM_PANELS * 2;

	// Default radius of the cylinder
	static const float DEFAULT_RADIUS;

private:
	// OpenGL vertex array object ID
	unsigned int vertexArrayID;

	// OpenGL vertex buffer object ID
	unsigned int vertexBufferID;

	// Array to store the vertex data for the cylinder
	float vertices[NUM_VERTICES];

	// Index to track the current position in the vertex array
	size_t verticesIndex = 0;

	void addVertexPoint(const Vec3 &point, const Vec3 &normal, bool isTop);
	void genCylinderVertices();

public:
	ConnectorTemplate();
	const float *getVerticesPtr() const;
	size_t getVerticesLength() const;
	void bind() const;
};
