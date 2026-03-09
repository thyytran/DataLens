#include "graphics/Window.h"

/*
 *  Constructor for the Window class.
 *  This constructor creates a GLFW window with the specified title, width, height, and resizability.
 *  @param windowTitle The title of the window.
 *  @param screenWidth The width of the window in pixels.
 *  @param screenHeight The height of the window in pixels.
 *  @param resizable A boolean indicating whether the window is resizable.
 */
Window::Window(const std::string &windowTitle, int screenWidth, int screenHeight, bool resizable) {
	if (!resizable) {
		glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
	}

	window = glfwCreateWindow(screenWidth, screenHeight, windowTitle.c_str(), 0, 0);
	if (!window) {
		std::cout << "Err > Window failed to initialize.\n\n";
		glfwTerminate();
	}
	glfwMakeContextCurrent(window);

	if (resizable) {
		glfwSetFramebufferSizeCallback(window, frameBufferSizeCallback);
	}

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_SCISSOR_TEST);
	glViewport(0, 0, screenWidth, screenHeight);
	glScissor(0, 0, screenWidth, screenHeight);
}

/*
 *  Framebuffer size callback function.
 *  This function is called when the framebuffer size of the window changes.
 *  It updates the viewport to match the new framebuffer size.
 *  @param window A pointer to the GLFW window.
 *  @param width The new width of the framebuffer in pixels.
 *  @param height The new height of the framebuffer in pixels.
 */
void Window::frameBufferSizeCallback(GLFWwindow *window, int width, int height) {
	glViewport(0, 0, width, height);
}

/*
 *  Sets the window's "should close" flag.
 *  @param shouldClose A boolean indicating whether the window should close.
 */
void Window::setShouldClose(bool shouldClose) {
	glfwSetWindowShouldClose(window, shouldClose);
}

/*
 *  Gets the window's "should close" flag.
 *  @return A boolean indicating whether the window should close.
 */
bool Window::shouldClose() const {
	return glfwWindowShouldClose(window);
}

/*
 *  Gets the GLFW window pointer.
 *  @return A pointer to the GLFW window.
 */
GLFWwindow *Window::get() const {
	return window;
}

/*
 *  Gets the width of the window in pixels.
 *  @return The width of the window in pixels.
 */
int Window::getWidth() const {
	int viewportData[4];
	glGetIntegerv(GL_VIEWPORT, viewportData);
	return viewportData[2];
}

/*
 *  Gets the height of the window in pixels.
 *  @return The height of the window in pixels.
 */
int Window::getHeight() const {
	int viewportData[4];
	glGetIntegerv(GL_VIEWPORT, viewportData);
	return viewportData[3];
}

/*
 *  Swaps the front and back buffers of the window.
 *  This function presents the rendered image to the screen.
 */
void Window::swapBuffers() {
	glfwSwapBuffers(window);
}

/*
 *  Clears the window with the specified color.
 *  @param r The red component of the color (0-255).
 *  @param g The green component of the color (0-255).
 *  @param b The blue component of the color (0-255).
 */
void Window::clear(int r, int g, int b) {
	const float MAX_COLOR_VALUE = 255.0f;
	glClearColor(
		r / MAX_COLOR_VALUE, 
		g / MAX_COLOR_VALUE, 
		b / MAX_COLOR_VALUE, 
		1.0f
	);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
