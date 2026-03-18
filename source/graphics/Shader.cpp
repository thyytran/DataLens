#include "graphics/Shader.h"

// Initialize static pointers to default shaders (sphere and connector) to nullptr
Shader *Shader::sphere = nullptr;
Shader *Shader::connector = nullptr;

/*
 *  Constructor for the Shader class.
 *  This constructor creates a shader program from vertex and fragment shader source code.
 *  @param vertexShaderSource   The source code of the vertex shader.
 *  @param fragmentShaderSource The source code of the fragment shader.
 */
Shader::Shader(const std::string &vertexShaderSource, const std::string &fragmentShaderSource) {
	readShaders(vertexShaderSource, fragmentShaderSource);
	projectionLocation = getUniformLocation("projection");
	viewLocation = getUniformLocation("view");
	modelLocation = getUniformLocation("model");
	normalLocation = getUniformLocation("normalMat");
}

/*
 * Constructor for the Shader class.
 * This constructor creates a shader program from source. 
 */
Shader::Shader(const std::string& vertexSource, const std::string& fragmentSource, bool fromSource) {
	id = 0;
	projectionLocation = 0;
	viewLocation = 0;
	modelLocation = 0;
	normalLocation = 0;

	// Compile vertex shader
	unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
	const char* vSource = vertexSource.c_str();
	glShaderSource(vertexShader, 1, &vSource, NULL);
	glCompileShader(vertexShader);
	checkShaderCompileErrors(vertexShader, ShaderType::VERTEX);

	// Compile fragment shader
	unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	const char* fSource = fragmentSource.c_str();
	glShaderSource(fragmentShader, 1, &fSource, NULL);
	glCompileShader(fragmentShader);
	checkShaderCompileErrors(fragmentShader, ShaderType::FRAGMENT);

	// Link program
	id = glCreateProgram();
	glAttachShader(id, vertexShader);
	glAttachShader(id, fragmentShader);
	glLinkProgram(id);
	checkShaderLinkErrors(id);

	// Clean up
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	// Get uniform locations
	projectionLocation = glGetUniformLocation(id, "u_projection");
	viewLocation = glGetUniformLocation(id, "u_view");
	modelLocation = glGetUniformLocation(id, "u_model");
	normalLocation = glGetUniformLocation(id, "u_normal");
}

/*
 *  Checks for shader compilation errors.
 *  @param id         The ID of the shader object.
 *  @param shaderType The type of the shader (vertex or fragment).
 */
void Shader::checkShaderCompileErrors(unsigned int id, ShaderType shaderType) {
	const int Err_BUFFER_SIZE = 512;
	int success;
	char infoLog[Err_BUFFER_SIZE];
	glGetShaderiv(id, GL_COMPILE_STATUS, &success);

	if (!success) {
		glGetShaderInfoLog(id, Err_BUFFER_SIZE, 0, infoLog);
		const char *shaderTypeStr = "Shader";
		if (shaderType == ShaderType::VERTEX) {
			shaderTypeStr = "Vertex shader";
		}
		else if (shaderType == ShaderType::FRAGMENT) {
			shaderTypeStr = "Fragment shader";
		}
		std::cerr << "Err > " << shaderTypeStr << " failed to compile: " << infoLog << "\n\n";
	}
}

/*
 *  Checks for shader linking errors.
 *  @param id The ID of the shader program.
 */
void Shader::checkShaderLinkErrors(unsigned int id) {
	int success;
	glGetProgramiv(id, GL_LINK_STATUS, &success);

	if (!success) {
		std::cerr << "Err > Shaders failed to link.\n\n";
	}
}

/*
 *  Reads and compiles the vertex and fragment shaders.
 *  @param vertexShader   The source code of the vertex shader.
 *  @param fragmentShader The source code of the fragment shader.
 */
void Shader::readShaders(const std::string &vertexShader, const std::string &fragmentShader) {
	const char* vertexShader_cstr = vertexShader.c_str();
	const char* fragmentShader_cstr = fragmentShader.c_str();

	//Compile Vertex Shader
	unsigned int vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShaderID, 1, &vertexShader_cstr, NULL);
	glCompileShader(vertexShaderID);
	checkShaderCompileErrors(vertexShaderID, ShaderType::VERTEX);

	//Compile Fragment Shader
	unsigned int fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShaderID, 1, &fragmentShader_cstr, NULL);
	glCompileShader(fragmentShaderID);
	checkShaderCompileErrors(fragmentShaderID, ShaderType::FRAGMENT);

	//Link Shaders
	unsigned int shaderProgramID = glCreateProgram();
	glAttachShader(shaderProgramID, vertexShaderID);
	glAttachShader(shaderProgramID, fragmentShaderID);
	glLinkProgram(shaderProgramID);
	this->id = shaderProgramID;
	checkShaderLinkErrors(id);

	glDeleteShader(vertexShaderID);
	glDeleteShader(fragmentShaderID);
}

/*
 *  Uses the shader program.
 */
void Shader::useProgram() const {
	glUseProgram(id);
}

/*
 *  Gets the location of a uniform variable in the shader program.
 *  @param uniformName The name of the uniform variable.
 *  @return The location of the uniform variable.
 */
unsigned int Shader::getUniformLocation(const std::string &uniformName) const {
	return glGetUniformLocation(id, uniformName.c_str());
}

/*
 *  Gets the ID of the shader program.
 *  @return The ID of the shader program.
 */
unsigned int Shader::getID() const {
	return id;
}

/*
 *  Sets the projection matrix uniform variable in the shader program.
 *  @param mat The projection matrix.
 */
void Shader::setProjectionMatrix(const Mat4 &mat) const {
	glUniformMatrix4fv(projectionLocation, 1, GL_FALSE, mat.getPtr());
}

/*
 *  Sets the view matrix uniform variable in the shader program.
 *  @param mat The view matrix.
 */
void Shader::setViewMatrix(const Mat4 &mat) const {
	glUniformMatrix4fv(viewLocation, 1, GL_FALSE, mat.getPtr());
}

/*
 *  Sets the model matrix uniform variable in the shader program.
 *  @param mat The model matrix.
 */
void Shader::setModelMatrix(const Mat4 &mat) const {
	glUniformMatrix4fv(modelLocation, 1, GL_FALSE, mat.getPtr());
}

/*
 *  Sets the normal matrix uniform variable in the shader program.
 *  @param mat The normal matrix.
 */
void Shader::setNormalMatrix(const Mat3 &mat) const {
	glUniformMatrix3fv(normalLocation, 1, GL_FALSE, mat.getPtr());
}

/*
 *  Sets a vec3 uniform variable in the shader program.
 *  @param uniformName The name of the uniform variable.
 *  @param value1      The first value of the vec3.
 *  @param value2      The second value of the vec3.
 *  @param value3      The third value of the vec3.
 */
void Shader::setVec3(
	const std::string &uniformName, 
	float value1, float value2, float value3
) const {
	unsigned int uniformLocation = getUniformLocation(uniformName);
	glUniform3f(uniformLocation, value1, value2, value3);
}

/*
 *  Loads the default shaders (sphere and connector).
 */
void Shader::loadDefaultShaders() {
	std::string sphereVertexShader = 
R"(#version 330 core

layout(location = 0) in vec3 vertexPos;
layout(location = 1) in vec3 sphereCenter;
layout(location = 2) in float sphereRadius;
layout(location = 3) in vec3 color;

out vec3 fragmentPos;
out vec3 fragmentNormal;
out vec3 fragmentColor;

uniform mat4 model;
uniform mat3 normalMat;
uniform mat4 view;
uniform mat4 projection;

void main() {
	vec4 worldPos = model * vec4(vertexPos * sphereRadius + sphereCenter, 1.0f);
	gl_Position = projection * view * worldPos;

	fragmentPos = vec3(worldPos);
	fragmentNormal = normalMat * vertexPos;
	fragmentColor = color;
}
)";
	std::string sphereFragmentShader = 
R"(#version 330 core

in vec3 fragmentPos;
in vec3 fragmentNormal;
in vec3 fragmentColor;

out vec4 finalFragmentColor;

void main() {
	vec3 lightColor = vec3(1.0f, 1.0f, 1.0f);
	vec3 lightDirection = vec3(0.0f, 0.0f, 1.0f);

	//Ambient
	float ambientStrength = 0.5f;
	vec3 ambient = ambientStrength * lightColor;

	//Diffuse
	vec3 normal = normalize(fragmentNormal);
	float diffuseStrength = max(dot(normal, lightDirection), 0.1f);
	vec3 diffuse = 0.65f * diffuseStrength * lightColor;

	finalFragmentColor = vec4((ambient + diffuse) * fragmentColor, 1.0f);
}
)";
	sphere = new Shader(sphereVertexShader, sphereFragmentShader);

	std::string connectorVertexShader =
R"(#version 330 core

layout(location = 0) in vec3 vertexPos;
layout(location = 1) in vec3 vertexNormal;
layout(location = 2) in float posIndicator;

out vec3 fragmentPos;
out vec3 fragmentNormal;
out vec3 fragmentColor;

uniform mat4 model;
uniform mat3 normalMat;
uniform mat4 view;
uniform mat4 projection;

uniform vec3 topColor;
uniform vec3 bottomColor;

void main() {
	vec4 worldPos = model * vec4(vertexPos, 1.0f);
	gl_Position = projection * view * worldPos;

	fragmentPos = vec3(worldPos);
	fragmentNormal = normalMat * vertexNormal;
	fragmentColor = bottomColor * posIndicator + topColor * (1.0f - posIndicator);
}
)";
	std::string connectorFragmentShader =
R"(#version 330 core

in vec3 fragmentPos;
in vec3 fragmentNormal;
in vec3 fragmentColor;

out vec4 finalFragmentColor;

void main() {
	vec3 lightColor = vec3(1.0f, 1.0f, 1.0f);
	vec3 lightDirection = vec3(0.0f, 0.0f, 1.0f);

	//Ambient
	float ambientStrength = 0.5f;
	vec3 ambient = ambientStrength * lightColor;

	//Diffuse
	vec3 normal = normalize(fragmentNormal);
	float diffuseStrength = max(dot(normal, lightDirection), 0.1f);
	vec3 diffuse = 0.65f * diffuseStrength * lightColor;

	finalFragmentColor = vec4((ambient + diffuse) * fragmentColor, 1.0f);
}
)";
	connector = new Shader(connectorVertexShader, connectorFragmentShader);
}

/*
 *  Gets the default sphere shader.
 *  @return A pointer to the default sphere shader.
 */
const Shader *Shader::getSphereDefault() {
	return sphere;
}

/*
 *  Gets the default connector shader.
 *  @return A pointer to the default connector shader.
 */
const Shader *Shader::getConnectorDefault() {
	return connector;
}

Shader* Shader::createFromSource(const char* vertexSource, const char* fragmentSource) {
	Shader* shader = new Shader();

	// Compile vertex shader
	unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShader, 1, &vertexSource, NULL);
	glCompileShader(vertexShader);
	shader->checkShaderCompileErrors(vertexShader, ShaderType::VERTEX);

	// Compile fragment shader
	unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
	glCompileShader(fragmentShader);
	shader->checkShaderCompileErrors(fragmentShader, ShaderType::FRAGMENT);

	// Link program
	shader->id = glCreateProgram();
	glAttachShader(shader->id, vertexShader);
	glAttachShader(shader->id, fragmentShader);
	glLinkProgram(shader->id);
	shader->checkShaderLinkErrors(shader->id);

	// Clean up
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);

	// Cache uniform locations
	shader->projectionLocation = glGetUniformLocation(shader->id, "u_projection");
	shader->viewLocation = glGetUniformLocation(shader->id, "u_view");
	shader->modelLocation = glGetUniformLocation(shader->id, "u_model");
	shader->normalLocation = glGetUniformLocation(shader->id, "u_normal");

	return shader;
}

/*
 *  Frees the resources used by the Shader class.
 */
void Shader::freeResources() {
	delete sphere;
	sphere = nullptr;

	delete connector;
	connector = nullptr;
}
