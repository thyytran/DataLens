#include "ResourceManager.h"

/*
 *  Initializes GLFW (Graphics Library Framework) with default OpenGL version.
 *  GLFW is used for window and input management.
 *  @return True if GLFW initialization was successful, false otherwise.
 */
bool ResourceManager::initGLFW() {
	if (!glfwInit()) {
		std::cout << "Err > GLFW failed to initialize.\n\n";
		return 0;
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, DEFAULT_MAJOR_GL);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, DEFAULT_MINOR_GL);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	return 1;
}

/*
 *  Initializes GLFW with a specified OpenGL version.
 *  This allows specifying a non-default OpenGL version for the application.
 *  @param majorOpenGLVersion The desired major OpenGL version.
 *  @param minorOpenGLVersion The desired minor OpenGL version.
 *  @return True if GLFW initialization was successful, false otherwise.
 */
bool ResourceManager::initGLFW(int majorOpenGLVersion, int minorOpenGLVersion) {
	if (!glfwInit()) {
		std::cout << "Err > GLFW failed to initialize.\n\n";
		return 0;
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, majorOpenGLVersion);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, minorOpenGLVersion);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	return 1;
}

/*
 *  Initializes GLEW (OpenGL Extension Wrangler Library).
 *  GLEW is used to access OpenGL extensions and functions.
 *  @return True if GLEW initialization was successful, false otherwise.
 */
bool ResourceManager::initGLEW() {
	GLenum error = glewInit();
	if (error != GLEW_OK) {
		std::cout << "Err > GLEW failed to initialize.\n\n";
		glfwTerminate();
		return 0;
	}
	return 1;
}

/*
 *  Initializes OpenGL settings.
 *  This sets up basic OpenGL states such as blending and face culling.
 */
void ResourceManager::initOpenGL() {
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_CULL_FACE);

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
}

/*
 *  Frees resources used by the ResourceManager.
 *  This terminates GLFW, releasing any resources it holds.
 */
void ResourceManager::freeResources() {
	glfwTerminate();
}
