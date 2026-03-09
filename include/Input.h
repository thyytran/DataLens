#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "graphics/Window.h"

// Enumeration for representing different keys on the keyboard
enum class Key {

	// Represents an undefined or invalid key
	UNDEFINED,

	// Lock Keys
	CAPS_LOCK, NUM_LOCK, SCROLL_LOCK,

	// Navigation Keys
	PAGE_UP, PAGE_DOWN, HOME, END, PRINT_SCREEN, PAUSE,

	// Control Keys
	LEFT_CTRL, RIGHT_CTRL, LEFT_ALT, RIGHT_ALT,

	// Special Keys
	SPACE, TAB, ENTER, BACKSPACE, DELETE_KEY, INSERT, ESCAPE, //DELETE is a macro defined in winnt.h
	
	// Arrow Keys
	UP, DOWN, LEFT, RIGHT,
	
	// Alphabet Keys
	A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z,
	
	// Number Keys (Main Row)
	MAIN_0, MAIN_1, MAIN_2, MAIN_3, MAIN_4, MAIN_5, MAIN_6, MAIN_7, MAIN_8, MAIN_9,
	
	// Punctuation Keys
	COMMMA, PERIOD, FORWARD_SLASH, BACK_SLASH, SEMICOLON, APOSTROPHE, LEFT_BRACKET, RIGHT_BRACKET, MINUS, EQUALS, BACKTICK,
	
	// Number Pad Keys
	PAD_0, PAD_1, PAD_2, PAD_3, PAD_4, PAD_5, PAD_6, PAD_7, PAD_8, PAD_9,
	
	// Number Pad Operators
	PAD_ADD, PAD_SUBTRACT, PAD_MULTIPLY, PAD_DIVIDE, PAD_DECIMAL,
	
	// Function Keys
	F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21, F22, F23, F24, F25
};

// Enumeration for representing different mouse buttons
enum class MouseButton {
	LEFT, MIDDLE, RIGHT
};

// Enumeration for representing different cursor types
enum class CursorType {
	UNDEFINED, NORMAL, HIDDEN, DISABLED
};

class Input {
private:
	// Number of keys being tracked
	const static unsigned int NUM_KEYS = 112;

	// Array mapping Key enum values to GLFW key codes
	const static int GLFW_KEY_LIST[NUM_KEYS];

	// Number of mouse buttons being tracked
	const static unsigned int NUM_MOUSE_BUTTONS = 3;

	// Array mapping MouseButton enum values to GLFW mouse button codes
	const static int GLFW_MOUSE_BUTTON_LIST[NUM_MOUSE_BUTTONS];

	// Static variable to store the current cursor type
	static CursorType cursorType;

	// Static variables to store the current mouse position
	static double mouseX, mouseY;

	// Static variables to store the mouse scroll offset
	static double mouseScrollOffsetX, mouseScrollOffsetY;

	static void mouseCursorCallback(GLFWwindow *window, double xPos, double yPos);
	static void mouseScrollCallback(GLFWwindow *window, double xOffset, double yOffset);

public:
	static bool keyPressed(Window *window, Key key);
	static bool mouseButtonPressed(Window *window, MouseButton mouseButton);
	static void initCallbacks(Window* window); //Just added
	static void pollInput(Window *window);

	static void setCursorType(Window *window, CursorType cursorType);

	static CursorType getCursorType();
	static double getMouseX();
	static double getMouseY();
	static double getMouseScrollOffsetX();
	static double getMouseScrollOffsetY();
};
