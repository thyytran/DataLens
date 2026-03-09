#pragma once

#include "graphics/Shader.h"
#include "graphics/Window.h"
#include "math/Vec.h"
#include "math/Mat.h"
#include "math/MathUtils.h"
// #include "../lib/glm/detail/type_vec.hpp"

const float YAW = -90.0f;
const float PITCH = 0.0f;
const float SPEED = 2.0f;
const float SENSITIVITY = 0.05f;
const float ZOOM = 45.0f;

class Camera {
private:
	const Vec3 START_POS;
	const float START_FOV = 60.0f;

	Vec3 position;
	float fov = 30.0f;

	// Just added
	Vec3 target;  // What camera looks at


	//void updateCameraVectors();

public:
	Camera(const Vec3 &initialPosition);
	void move(const Vec3 &values);
	void zoom(float value);
	void reset();
	void Reset(Vec3 position, float yaw = -90, float pitch = 0, float targetDistance = 5);
	Mat4 getViewMatrix() const;
	Mat4 getProjectionMatrix(const Window *window) const;

	// Just added
	const Vec3& getPosition() const { return position; }
	float getFOV() const { return fov; }
	void setTarget(const Vec3& t) { target = t; }

	// Euler angles
	float Yaw;
	float Pitch;
	float Theta = 90;
	float Phi{};

	// Camera options
	float MovementSpeed;
	float MouseSensitivity;
	float Zoom;
	float ZoomSmooth{};

	// Camera attributes
	// glm::vec3 Position, OrbitPosition;
	// glm::vec3 Forward = glm::vec3(0.0f, 0.0f, -1.0f);
	// glm::vec3 Up;
	// glm::vec3 Right;
	// glm::vec3 WorldUp = glm::vec3(0.0f, 1.0f, 0.0f);

	// Camera attributes - using custom Vec libraries
	Vec3 Position, OrbitPosition;
	Vec3 Forward = Vec3(0.0f, 0.0f, -1.0f);
	Vec3 Up;
	Vec3 Right;
	Vec3 WorldUp = Vec3(0.0f, 1.0f, 0.0f);
	
	Vec3 _position;

	// glm::vec3 Target;
	// glm::vec3 TargetSmooth;
	float TargetDistance = 5.0f;
	float TargetDistanceSmooth = TargetDistance;

	Vec3 Target;
	Vec3 TargetSmooth;
};
