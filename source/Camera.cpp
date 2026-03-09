#include "Camera.h"
// #include "../lib/glm/detail/func_trigonometric.hpp"
// #include "../lib/glm/detail/func_geometric.hpp"

/*
void Camera::updateCameraVectors()
{
	Forward.setX(cos(glm::radians(Yaw)) * cos(glm::radians(Pitch)));
	Forward.setY(sin(glm::radians(Pitch)));
	Forward.setZ(sin(glm::radians(Yaw)) * cos(glm::radians(Pitch)));
	// Forward = glm::normalize(Forward); // 
	// Forward = MathUtils

	// also re-calculate the Right and Up vector
	//Right = glm::normalize(glm::cross(Forward, WorldUp));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
	//Up = glm::normalize(glm::cross(Right, Forward));
}
*/


/*
 *  Constructor for the Camera class.
 *  @param initialPosition The initial position of the camera in 3D space.
 */

//Camera::Camera(const Vec3 &initialPosition) :
//	START_POS(initialPosition), position(initialPosition) {}


Camera::Camera(const Vec3& initialPosition) :
	START_POS(initialPosition), position(initialPosition), target(0, 0, 0) {
}  // Add target

/*
 *  Moves the camera by adding the given values to its current position.
 *  @param values A Vec3 representing the amount to move the camera in each axis.
 */ 
void Camera::move(const Vec3 &values) {
	position += values;
}

/*
 *  Zooms the camera by changing its field of view (FOV).
 *  @param value The amount to change the FOV. Positive values zoom out, negative values zoom in.
 */
void Camera::zoom(float value) {
	fov -= value;
	if (fov < 0.1f) {
		fov = 0.1f;
	}
	else if (fov > 179.0f) {
		fov = 179.0f;
	}
}

/*
 *  Resets the camera to its initial position and FOV.
 */
void Camera::reset() {
	position = START_POS;
	fov = START_FOV;
}

/*
 *  Gets the view matrix for the camera.
 *  The view matrix transforms world coordinates into camera coordinates.
 *  @return A Mat4 representing the view matrix.
 */
/*
Mat4 Camera::getViewMatrix() const {
	return MathUtils::MatGen::lookAt(
		position,
//		Vec3(position.getX(), position.getY(), 0.0f), // 1.0f to 0.0f
		Vec3(0.0f, 0.0f, 0.0f), // Look at origin (where the protein is)
		Vec3(0.0f, 1.0f, 0.0f)
	);
}
*/

Mat4 Camera::getViewMatrix() const {
	return MathUtils::MatGen::lookAt(
		position,
		target,  //  Look at target, not fixed origin
		Vec3(0.0f, 1.0f, 0.0f)
	);
}

/*
 *  Gets the projection matrix for the camera.
 *  The projection matrix transforms camera coordinates into clip space coordinates.
 *  @param window A pointer to the Window object, used to get the window dimensions.
 *  @return A Mat4 representing the projection matrix.
 */
Mat4 Camera::getProjectionMatrix(const Window *window) const {
	return MathUtils::MatGen::perspective<float>(
		MathUtils::toRadians(fov),
		(float)window->getWidth(),
		(float)window->getHeight(),
		2.0f,
		20.0f
	);
}

/*
 *  Resets the camera to a specified position, yaw, pitch, and target distance.
 *  @param position       The new position of the camera.
 *  @param yaw            The new yaw (horizontal angle) of the camera.
 *  @param pitch          The new pitch (vertical angle) of the camera.
 *  @param targetDistance The distance from the camera to the target.
 */
void Camera::Reset(Vec3 position, float yaw, float pitch, float targetDistance)
{
	Yaw = yaw;
	Pitch = pitch;
	Theta = pitch + 90;
	Phi = 0;
	Zoom = 45;
	//Target = position;
	//TargetDistance = targetDistance;
	//Position = position - TargetDistance * Forward;
	//updateCameraVectors();
}

