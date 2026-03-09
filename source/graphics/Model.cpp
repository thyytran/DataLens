#include "graphics/Model.h"

unsigned int Model::vertexArrayID;
unsigned int Model::vertexBufferID;
unsigned int Model::sphereBufferID;
Mat4 Model::modelMatrix(MathUtils::MatGen::identity<float, 4>());

std::vector<float> Model::atomSpheres;
std::vector<float> Model::backboneSpheres;
std::vector<float> Model::disulfideBondSpheres;

std::vector<Connector*> Model::backboneSegments;
std::vector<Connector*> Model::disulfideBonds;

const SphereTemplate *Model::sphereTemplate = nullptr;
const ConnectorTemplate *Model::connectorTemplate = nullptr;

const MoleculeData *Model::moleculeData = nullptr;

/*
 *  Sets the attributes for the sphere buffer.
 *  This function defines how the data in the sphere buffer is interpreted by the shader.
 */
void Model::setSphereBufferAttributes() {
	//Center
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (GLvoid*)0);
	glVertexAttribDivisor(1, 1);
	glEnableVertexAttribArray(1);

	//Radius
	glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
	glVertexAttribDivisor(2, 1);
	glEnableVertexAttribArray(2);

	//Color
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(GLfloat), (GLvoid*)(4 * sizeof(GLfloat)));
	glVertexAttribDivisor(3, 1);
	glEnableVertexAttribArray(3);
}

/*
 *  Adds a sphere to the specified vector.
 *  @param center The center coordinates of the sphere.
 *  @param radius The radius of the sphere.
 *  @param color  A pointer to the color of the sphere.
 *  @param vec    A pointer to the vector to add the sphere data to.
 */
void Model::addSphere(
	const Vec3 &center,
	float radius,
	const Color *color,
	std::vector<float> *vec
) {
	vec->reserve(vec->size() + 7);
	vec->push_back(center.getX());
	vec->push_back(center.getY());
	vec->push_back(center.getZ());
	vec->push_back(radius);
	vec->push_back(color->r);
	vec->push_back(color->g);
	vec->push_back(color->b);
}

/*
 *  Checks if the selection is valid.
 *  @param selection A pointer to the selection object.
 *  @return True if the selection is valid, false otherwise.
 */
bool Model::selectionIsValid(const Selection *selection) {
	if (!moleculeData || moleculeData->atoms.size() == 0) {
		std::cerr << "Err > No atoms found\n\n";
		return false;
	}
	return true;
}

/*
 *  Colors the atoms based on the provided function and selection.
 *  @param func      A function that takes an Atom pointer as input and returns a Color.
 *  @param selection A pointer to the selection object.
 */
void Model::colorAtomsByFunc(
	std::function<Color(const Atom*)> func,
	const Selection *selection
) {
	if (!selectionIsValid(selection)) {
		return;
	}
	for (size_t i = 4; i <= atomSpheres.size() - 3; i += 7) {
		const Atom *atom = &moleculeData->atoms[(i - 4) / 7];
		if (selection->isMatch(atom)) {
			Color color;
			try {
				color = func(atom);
			}
			catch (const std::invalid_argument &e) {
				//For invalid element names
				std::cerr << e.what() << "\n\n";
				continue;
			}
			atomSpheres[i] = color.r;
			atomSpheres[i + 1] = color.g;
			atomSpheres[i + 2] = color.b;
		}
	}
	syncSphereBuffer(true, false, false);
}

/*
 *  Colors the connectors based on the provided function, selection, and connector type.
 *  @param func          A function that takes a Connector pointer as input and returns a pair of Colors (top and bottom).
 *  @param selection     A pointer to the selection object.
 *  @param connectorType The type of connector to color (backbone or disulfide bond).
 */
void Model::colorConnectorsByFunc(
	std::function<std::pair<Color, Color>(const Connector*)> func,
	const Selection *selection,
	ConnectorType connectorType
) {
	std::vector<Connector*> *connectors = nullptr;
	std::vector<float> *spheres = nullptr;

	bool syncBackboneSpheres = false;
	bool syncDisulfideBondSpheres = false;

	if (connectorType == ConnectorType::BACKBONE) {
		connectors = &backboneSegments;
		spheres = &backboneSpheres;
		syncBackboneSpheres = true;
	}
	else {
		connectors = &disulfideBonds;
		spheres = &disulfideBondSpheres;
		syncDisulfideBondSpheres = true;
	}

	for (size_t i = 0; i < connectors->size(); ++i) {
		bool atom1IsMatch = selection->isMatch((*connectors)[i]->atom1);
		bool atom2IsMatch = selection->isMatch((*connectors)[i]->atom2);
		if (atom1IsMatch || atom2IsMatch) {
			std::pair<Color, Color> colors = func((*connectors)[i]);
			try {
				colors = func((*connectors)[i]);
			}
			catch (const std::invalid_argument &e) {
				//For invalid element names
				std::cerr << e.what() << "\n\n";
				continue;
			}

			if (atom1IsMatch) {
				(*connectors)[i]->setColors(&colors.first, nullptr);

				//Set color of top connector sphere
				size_t sphereIndex = i * 2;
				size_t sphereColorIndex = sphereIndex * 7 + 4;

				(*spheres)[sphereColorIndex] = colors.first.r;
				(*spheres)[sphereColorIndex + 1] = colors.first.g;
				(*spheres)[sphereColorIndex + 2] = colors.first.b;
			}
			if (atom2IsMatch) {
				(*connectors)[i]->setColors(nullptr, &colors.second);

				//Set color of top connector sphere
				size_t sphereIndex = i * 2 + 1;
				size_t sphereColorIndex = sphereIndex * 7 + 4;

				(*spheres)[sphereColorIndex] = colors.second.r;
				(*spheres)[sphereColorIndex + 1] = colors.second.g;
				(*spheres)[sphereColorIndex + 2] = colors.second.b;
			}
		}
	}
	syncSphereBuffer(false, syncBackboneSpheres, syncDisulfideBondSpheres);
}

/*
 *  Initializes auto-rotation parameters (currently empty).
 */
void Model::initAutoRotation()
{
	//rotationSpeed = { 0.0f, 0.0f, 0.0f };
	//rotationAccumulator = 0.0f;
}

/*
 *  Resets the model data.
 *  This function clears all the data associated with the model, including atom spheres,
 *  backbone segments, disulfide bonds, and molecule data.
 */
void Model::reset() {
	atomSpheres.clear();
	backboneSpheres.clear();
	disulfideBondSpheres.clear();

	for (size_t i = 0; i < backboneSegments.size(); ++i) {
		delete backboneSegments[i];
	}
	backboneSegments.clear();

	for (size_t i = 0; i < disulfideBonds.size(); ++i) {
		delete disulfideBonds[i];
	}
	disulfideBonds.clear();

	moleculeData = nullptr;
}

/*
 *  Sets the sphere template.
 *  @param sphereTemplate A pointer to the sphere template object.
 */
void Model::setSphereTemplate(const SphereTemplate *sphereTemplate) {
	Model::sphereTemplate = sphereTemplate;
}

/*
 *  Sets the connector template.
 *  @param connectorTemplate A pointer to the connector template object.
 */
void Model::setConnectorTemplate(const ConnectorTemplate *connectorTemplate) {
	Model::connectorTemplate = connectorTemplate;
}

/*
 *  Loads the molecule data.
 *  This function loads the molecule data, calculates the scaling and translation,
 *  creates the atom spheres and connectors, and fills the sphere buffer.
 *  @param moleculeData A pointer to the molecule data object.
 */
void Model::loadMoleculeData(const MoleculeData *moleculeData) {
	reset();

	//Gather coordinate information to adjust model size and position
	Vec3 coordSums;
	float maxCoord = 0.0f;
	for (size_t i = 0; i < moleculeData->atoms.size(); ++i) {
		coordSums += moleculeData->atoms[i].coords;
		for (unsigned int j = 0; j < 3; j++) {
			float absCoord = std::abs(moleculeData->atoms[i].coords.get(j));
			if (absCoord > maxCoord) {
				maxCoord = absCoord;
			}
		}
	}
	float coordScale = 7.5f / maxCoord;
	Vec3 coordAverages = coordSums * coordScale / (float)moleculeData->atoms.size();

	std::vector<const Atom*> alphaCarbons;

	for (size_t i = 0; i < moleculeData->atoms.size(); ++i) {
		//Add to alpha carbon data for constructing backbone
		if (moleculeData->atoms[i].name == "CA") {
			alphaCarbons.push_back(&(moleculeData->atoms[i]));
		}

		Color color;
		try {
			color = Color::fromElement(moleculeData->atoms[i].element);
		}
		catch (const std::invalid_argument &e) {
			std::cout << e.what() << "\n\n";
			continue;
		}
		
		addAtom(
			moleculeData->atoms[i].coords * coordScale - coordAverages,
			SphereTemplate::DEFAULT_RADIUS,
			&color
		);
	}

	//Construct backbone using connectors
	if (alphaCarbons.size() >= 2) {
		for (size_t i = 0; i < alphaCarbons.size() - 1; ++i) {
			//Same chain, no residues skipped
			if (alphaCarbons[i]->chain == alphaCarbons[i + 1]->chain &&
				(alphaCarbons[i]->residueNum == alphaCarbons[i + 1]->residueNum ||
					alphaCarbons[i]->residueNum + 1 == alphaCarbons[i + 1]->residueNum)) {

				Vec3 coords1 = alphaCarbons[i]->coords * coordScale - coordAverages;
				Vec3 coords2 = alphaCarbons[i + 1]->coords * coordScale - coordAverages;
				Color color = Color::fromElement("C");

				addConnector(
					alphaCarbons[i], alphaCarbons[i + 1],
					ConnectorTemplate::DEFAULT_RADIUS,
					&color, &color, coords1, coords2,
					ConnectorType::BACKBONE
				);
			}
		}
	}

	//Add disulfide bonds
	for (size_t i = 0; i < moleculeData->disulfideBonds.size(); ++i) {
		//Start residue same as end residue
		if (moleculeData->disulfideBonds[i].chain1 == moleculeData->disulfideBonds[i].chain2 &&
			moleculeData->disulfideBonds[i].residue1 == moleculeData->disulfideBonds[i].residue2) {

			continue;
		}
		
		const Atom *atom1 = nullptr;
		const Atom *atom2 = nullptr;

		//Find sulfur atoms
		for (size_t j = 0; j < moleculeData->atoms.size(); ++j) {
			if (atom1 && atom2) {
				break;
			}
			
			if (!atom1 &&
				moleculeData->atoms[j].chain == moleculeData->disulfideBonds[i].chain1 &&
				moleculeData->atoms[j].residueNum == moleculeData->disulfideBonds[i].residue1 &&
				moleculeData->atoms[j].element == "S") {

				atom1 = &(moleculeData->atoms[j]);
			}

			if (!atom2 &&
				moleculeData->atoms[j].chain == moleculeData->disulfideBonds[i].chain2 &&
				moleculeData->atoms[j].residueNum == moleculeData->disulfideBonds[i].residue2 &&
				moleculeData->atoms[j].element == "S") {

				atom2 = &(moleculeData->atoms[j]);
			}
		}

		if (atom1 && atom2) {
			Vec3 coords1 = atom1->coords * coordScale - coordAverages;
			Vec3 coords2 = atom2->coords * coordScale - coordAverages;
			Color color = Color::fromElement("S");

			addConnector(
				atom1, atom2,
				0.0f, //Invisible by default
				&color, &color, coords1, coords2,
				ConnectorType::DISULFIDE_BOND
			);
		}
	}

	fillSphereBuffer();
	Model::moleculeData = moleculeData;
}

void Model::addAtom(const Vec3 &center, float radius, const Color *color) {
	addSphere(center, radius, color, &atomSpheres);
}

/*
 *  Adds a connector to the specified connector and sphere vectors.
 *  This function creates a new Connector object and adds it to the appropriate connector vector
 *  (backboneSegments or disulfideBonds) based on the connectorType. It also adds two spheres
 *  to the corresponding sphere vector (backboneSpheres or disulfideBondSpheres), representing
 *  the endpoints of the connector.
 *  @param atom1         A pointer to the first atom.
 *  @param atom2         A pointer to the second atom.
 *  @param radius        The radius of the connector.
 *  @param topColor      A pointer to the color of the top of the connector.
 *  @param bottomColor   A pointer to the color of the bottom of the connector.
 *  @param point1        The 3D coordinates of the first atom.
 *  @param point2        The 3D coordinates of the second atom.
 *  @param connectorType The type of connector (backbone or disulfide bond).
 */
void Model::addConnector(
	const Atom *atom1, const Atom *atom2,
	float radius,
	const Color *topColor, const Color *bottomColor,
	const Vec3 &point1, const Vec3 &point2,
	ConnectorType connectorType
) {
	std::vector<Connector*> *connectors = nullptr;
	std::vector<float> *spheres = nullptr;

	if (connectorType == ConnectorType::BACKBONE) {
		connectors = &backboneSegments;
		spheres = &backboneSpheres;
	}
	else {
		connectors = &disulfideBonds;
		spheres = &disulfideBondSpheres;
	}

	connectors->push_back(
		new Connector(
			atom1, atom2,
			radius,
			topColor, bottomColor,
			point1, point2
		)
	);

	addSphere(point1, radius, topColor, spheres);
	addSphere(point2, radius, bottomColor, spheres);
}

/*
 *  Generates the OpenGL buffers.
 *  This function generates the vertex array object (VAO) and the vertex buffer objects (VBOs)
 *  for the model.
 */
void Model::genBuffers() {
	glGenVertexArrays(1, &vertexArrayID);
	glGenBuffers(1, &vertexBufferID);
	glGenBuffers(1, &sphereBufferID);
}

/*
 *  Fills the sphere template buffer with data from the SphereTemplate object.
 *  This function binds the vertex array object and the vertex buffer object,
 *  and then copies the vertex data from the SphereTemplate object to the vertex buffer object.
 *  It also sets up the vertex attribute pointer for the vertex positions.
 */
void Model::fillSphereTemplateBuffer() {
	if (!sphereTemplate) {
		std::cerr << "Err > Sphere template not set\n\n";
		return;
	}

	const size_t TEMPLATE_SIZE = sizeof(float) * SphereTemplate::NUM_VERTICES;

	glBindVertexArray(vertexArrayID);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBufferID);

	glBufferData(GL_ARRAY_BUFFER, TEMPLATE_SIZE, sphereTemplate->getVerticesPtr(), GL_STATIC_DRAW);

	//Positions
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(0);
}

/*
 *  Allocates memory for the sphere buffer.
 *  This function calculates the total size of the sphere data (atom spheres, backbone spheres,
 *  and disulfide bond spheres) and allocates memory for the sphere buffer object.
 */
void Model::allocateSphereBuffer() {
	const size_t ATOM_SPHERES_SIZE = sizeof(float) * atomSpheres.size();
	const size_t BACKBONE_SPHERES_SIZE = sizeof(float) * backboneSpheres.size();
	const size_t DISULFIDE_BOND_SPHERES_SIZE = sizeof(float) * disulfideBondSpheres.size();
	const size_t SPHERES_SIZE = 
		ATOM_SPHERES_SIZE + BACKBONE_SPHERES_SIZE + DISULFIDE_BOND_SPHERES_SIZE;

	glBindVertexArray(vertexArrayID);
	glBindBuffer(GL_ARRAY_BUFFER, sphereBufferID);

	glBufferData(GL_ARRAY_BUFFER, SPHERES_SIZE, 0, GL_STATIC_DRAW);
	setSphereBufferAttributes();
}

/*
 *  Synchronizes the sphere buffer with the data in the sphere vectors.
 *  This function copies the data from the atomSpheres, backboneSpheres, and disulfideBondSpheres
 *  vectors to the sphere buffer object. It allows synchronizing only specific parts of the buffer
 *  by using the syncAtomSpheres, syncBackboneSpheres, and syncDisulfideBondSpheres flags.
 *  @param syncAtomSpheres       A flag indicating whether to synchronize the atom spheres data.
 *  @param syncBackboneSpheres   A flag indicating whether to synchronize the backbone spheres data.
 *  @param syncDisulfideBondSpheres A flag indicating whether to synchronize the disulfide bond spheres data.
 */
void Model::syncSphereBuffer(
	bool syncAtomSpheres, 
	bool syncBackboneSpheres,
	bool syncDisulfideBondSpheres
) {
	const size_t ATOM_SPHERES_SIZE = sizeof(float) * atomSpheres.size();
	const size_t BACKBONE_SPHERES_SIZE = sizeof(float) * backboneSpheres.size();
	const size_t DISULFIDE_BOND_SPHERES_SIZE = sizeof(float) * disulfideBondSpheres.size();

	glBindVertexArray(vertexArrayID);
	glBindBuffer(GL_ARRAY_BUFFER, sphereBufferID);

	if (syncAtomSpheres) {
		glBufferSubData(
			GL_ARRAY_BUFFER,
			0,
			ATOM_SPHERES_SIZE,
			atomSpheres.data()
		);
	}

	if (syncBackboneSpheres) {
		glBufferSubData(
			GL_ARRAY_BUFFER,
			ATOM_SPHERES_SIZE,
			BACKBONE_SPHERES_SIZE,
			backboneSpheres.data()
		);
	}

	if (syncDisulfideBondSpheres) {
		glBufferSubData(
			GL_ARRAY_BUFFER,
			ATOM_SPHERES_SIZE + BACKBONE_SPHERES_SIZE,
			DISULFIDE_BOND_SPHERES_SIZE,
			disulfideBondSpheres.data()
		);
	}
}

/*
 *  Fills the sphere buffer with data from the sphere vectors.
 *  This function copies the data from the atomSpheres, backboneSpheres, and disulfideBondSpheres
 *  vectors to the sphere buffer object.
 */
void Model::fillSphereBuffer() {
	const size_t ATOM_SPHERES_SIZE = sizeof(float) * atomSpheres.size();
	const size_t BACKBONE_SPHERES_SIZE = sizeof(float) * backboneSpheres.size();
	const size_t DISULFIDE_BOND_SPHERES_SIZE = sizeof(float) * disulfideBondSpheres.size();
	const size_t SPHERES_SIZE = 
		ATOM_SPHERES_SIZE + BACKBONE_SPHERES_SIZE + DISULFIDE_BOND_SPHERES_SIZE;

	glBindVertexArray(vertexArrayID);
	glBindBuffer(GL_ARRAY_BUFFER, sphereBufferID);

	glBufferData(GL_ARRAY_BUFFER, SPHERES_SIZE, 0, GL_STATIC_DRAW);
	setSphereBufferAttributes();

	//Atom spheres
	glBufferSubData(
		GL_ARRAY_BUFFER,
		0,
		ATOM_SPHERES_SIZE,
		atomSpheres.data()
	);

	//Backbone spheres
	glBufferSubData(
		GL_ARRAY_BUFFER,
		ATOM_SPHERES_SIZE,
		BACKBONE_SPHERES_SIZE,
		backboneSpheres.data()
	);

	//Disulfide bond spheres
	glBufferSubData(
		GL_ARRAY_BUFFER,
		ATOM_SPHERES_SIZE + BACKBONE_SPHERES_SIZE,
		DISULFIDE_BOND_SPHERES_SIZE,
		disulfideBondSpheres.data()
	);
}

/*
 *  Renders the model.
 *  This function renders the atom spheres and connectors using the provided shaders, window, and camera.
 *  @param sphereShader    A pointer to the shader program used for rendering the atom spheres.
 *  @param connectorShader A pointer to the shader program used for rendering the connectors.
 *  @param window          A pointer to the window object.
 *  @param camera          A pointer to the camera object.
 */
void Model::render(
	const Shader *sphereShader,
	const Shader *connectorShader,
	const Window *window,
	const Camera *camera
) {
	Mat4 viewMatrix = camera->getViewMatrix();
	Mat4 projectionMatrix = camera->getProjectionMatrix(window);

	//Draw atom and connector spheres
	if (!sphereShader) {
		std::cerr << "Err > Unloaded sphere shaders\n\n";
		return;
	}

	glBindVertexArray(vertexArrayID);
	sphereShader->useProgram();

	sphereShader->setModelMatrix(modelMatrix);
	sphereShader->setNormalMatrix(MathUtils::MatGen::normal<float>(modelMatrix));
	sphereShader->setViewMatrix(viewMatrix);
	sphereShader->setProjectionMatrix(projectionMatrix);

	glDrawArraysInstanced(
		GL_TRIANGLES,
		0,
		GLsizei(SphereTemplate::NUM_VERTICES / 3),
		GLsizei(
			(atomSpheres.size() + backboneSpheres.size() + disulfideBondSpheres.size()) / 7
		)
	);

	//Draw connectors
	if (!connectorShader) {
		std::cerr << "Err > Unloaded connector shaders\n\n";
		return;
	}
	if (!connectorTemplate) {
		std::cerr << "Err > Connector template not set\n\n";
		return;
	}

	connectorTemplate->bind();
	connectorShader->useProgram();

	connectorShader->setViewMatrix(viewMatrix);
	connectorShader->setProjectionMatrix(projectionMatrix);

	//Set connector model and normal matrices, set shader's color, draw arrays
	for(size_t i = 0; i < backboneSegments.size(); ++i) {
		backboneSegments[i]->render(connectorShader, modelMatrix, connectorTemplate);
	}
	for (size_t i = 0; i < disulfideBonds.size(); ++i) {
		disulfideBonds[i]->render(connectorShader, modelMatrix, connectorTemplate);
	}
}

/*
 *  Sets the radius of the atoms based on the provided selection.
 *  @param radius    The new radius of the atoms.
 *  @param selection A pointer to the selection object.
 *  @param reversed  A flag indicating whether to reverse the selection.
 */
void Model::setAtomRadius(float radius, const Selection *selection, bool reversed) {
	if (!selectionIsValid(selection)) {
		return;
	}
	for (size_t i = 3; i < atomSpheres.size(); i += 7) {
		if (selection->isMatch(&moleculeData->atoms[(i - 3) / 7], reversed)) {
			atomSpheres[i] = radius;
		}
	}
	syncSphereBuffer(true, false, false);
}

/*
 *  Sets the radius of the connectors based on the provided selection and connector type.
 *  @param radius        The new radius of the connectors.
 *  @param selection     A pointer to the selection object.
 *  @param connectorType The type of connector to set the radius for (backbone or disulfide bond).
 *  @param reversed      A flag indicating whether to reverse the selection.
 */
void Model::setConnectorRadius(
	float radius, const Selection *selection, ConnectorType connectorType,
	bool reversed
) {
	std::vector<Connector*> *connectors = nullptr;
	std::vector<float> *spheres = nullptr;
	
	bool syncBackboneSpheres = false;
	bool syncDisulfideBondSpheres = false;

	if (connectorType == ConnectorType::BACKBONE) {
		connectors = &backboneSegments;
		spheres = &backboneSpheres;
		syncBackboneSpheres = true;
	}
	else {
		connectors = &disulfideBonds;
		spheres = &disulfideBondSpheres;
		syncDisulfideBondSpheres = true;
	}

	for (size_t i = 0; i < connectors->size(); ++i) {
		bool atom1IsMatch = selection->isMatch((*connectors)[i]->atom1, reversed);
		bool atom2IsMatch = selection->isMatch((*connectors)[i]->atom2, reversed);
		if ((!reversed &&  atom1IsMatch && atom2IsMatch) ||
			(reversed && (atom1IsMatch || atom2IsMatch))) {

			(*connectors)[i]->setRadius(radius);

			//Set radius of connector spheres
			size_t sphereIndex1 = i * 2;
			size_t sphereIndex2 = i * 2 + 1;

			size_t sphereRadiusIndex1 = sphereIndex1 * 7 + 3;
			size_t sphereRadiusIndex2 = sphereIndex2 * 7 + 3;

			(*spheres)[sphereRadiusIndex1] = radius;
			(*spheres)[sphereRadiusIndex2] = radius;
		}
	}
	syncSphereBuffer(false, syncBackboneSpheres, syncDisulfideBondSpheres);
}

/*
 *  Sets the color of the atoms based on the provided color and selection.
 *  @param color     A pointer to the color to set for the atoms.
 *  @param selection A pointer to the selection object.
 */
void Model::setAtomColor(const Color *color, const Selection *selection) {
	colorAtomsByFunc(
		[&color](const Atom *atom) {
			return *color;
		},
		selection
	);
}

/*
 *  Colors the atoms based on their element.
 *  @param selection A pointer to the selection object.
 */
void Model::colorAtomsDefault(const Selection *selection) {
	colorAtomsByFunc(
		[](const Atom *atom) {
			return Color::fromElement(atom->element);
		},
		selection
	);
}

/*
 *  Colors the atoms based on their structure.
 *  @param selection A pointer to the selection object.
 */
void Model::colorAtomsByStructure(const Selection *selection) {
	colorAtomsByFunc(
		[&moleculeData = std::as_const(moleculeData)](const Atom *atom) {
			return Color::fromStructure(atom, moleculeData);
		},
		selection
	);
}

/*
 *  Sets the color of the connectors based on the provided color and selection.
 *  @param color         A pointer to the color to set for the connectors.
 *  @param selection     A pointer to the selection object.
 *  @param connectorType The type of connector to set the color for (backbone or disulfide bond).
 */
void Model::setConnectorColor(
	const Color *color, const Selection *selection,
	ConnectorType connectorType
) {
	colorConnectorsByFunc(
		[&color](const Connector *connector) {
			return std::make_pair(*color, *color);
		},
		selection, connectorType
	);
}

/*
 *  Colors the connectors based on their default color (carbon for backbone, sulfur for disulfide bond).
 *  @param selection     A pointer to the selection object.
 *  @param connectorType The type of connector to color (backbone or disulfide bond).
 */
void Model::colorConnectorsDefault(
	const Selection *selection, ConnectorType connectorType
) {
	colorConnectorsByFunc(
		[&connectorType](const Connector *connector) {
			if (connectorType == ConnectorType::BACKBONE) {
				Color color = Color::fromElement("C");
				return std::make_pair(color, color);
			}
			else {
				Color color = Color::fromElement("S");
				return std::make_pair(color, color);
			}
		},
		selection, connectorType
	);
}

/*
 *  Colors the connectors based on the structure of the atoms they connect.
 *  @param selection     A pointer to the selection object.
 *  @param connectorType The type of connector to color (backbone or disulfide bond).
 */
void Model::colorConnectorsByStructure(
	const Selection *selection, ConnectorType connectorType
) {
	colorConnectorsByFunc(
		[&moleculeData = std::as_const(moleculeData)](const Connector *connector) {
			Color color1 = Color::fromStructure(connector->atom1, moleculeData);
			Color color2 = Color::fromStructure(connector->atom2, moleculeData);
			return std::make_pair(color1, color2);
		},
		selection, connectorType
	);
}

/*
 *  Rotates the model by the specified angle in radians.
 *  @param angleRadians A Vec3 object representing the rotation angle in radians around each axis (x, y, z).
 */
void Model::rotate(const Vec3 &angleRadians) {
	modelMatrix = modelMatrix * MathUtils::MatGen::rotation<float>(angleRadians);
}

