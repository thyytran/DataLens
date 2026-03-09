#include "graphics/RayCaster.h"
#include "math/MathUtils.h"
#include <cmath>


Ray RayCaster::screenToRay(
    double mouseX, double mouseY,
    const Window* window,
    const Camera* camera
) {
    float screenWidth = (float)window->getWidth();
    float screenHeight = (float)window->getHeight();

    // Get camera position as ray origin
    Vec3 rayOrigin = camera->getPosition();

    // Convert screen coords to NDC [-1, 1]
    float x = (2.0f * (float)mouseX) / screenWidth - 1.0f;
    float y = 1.0f - (2.0f * (float)mouseY) / screenHeight;  // Flip Y

    // Get camera FOV and aspect ratio
    float fov = camera->getFOV();
    float aspect = screenWidth / screenHeight;

    // Calculate ray direction in view space
    float tanHalfFov = std::tan(fov * 0.5f * 3.14159f / 180.0f);
    float viewX = x * aspect * tanHalfFov;
    float viewY = y * tanHalfFov;
    float viewZ = -1.0f;  // Looking down -Z in view space

    Vec3 rayDirView(viewX, viewY, viewZ);

    // Transform ray direction to world space
    // Since your camera just orbits around origin, we can simplify:
    // The "forward" direction is from camera position toward origin
    Vec3 forward = (Vec3(0.0f, 0.0f, 0.0f) - rayOrigin).normalize();
    Vec3 right = Vec3(0.0f, 1.0f, 0.0f).cross(forward).normalize();
    Vec3 up = forward.cross(right);

    // Build ray direction in world space
    Vec3 rayDirection = (
        right * viewX +
        up * viewY +
        forward
        ).normalize();

    return Ray(rayOrigin, rayDirection);
}

std::optional<float> RayCaster::raySphereIntersection(
    const Ray& ray,
    const Vec3& sphereCenter,
    float sphereRadius
) {
    Vec3 oc = ray.origin - sphereCenter;

    float a = ray.direction.dot(ray.direction);
    float b = 2.0f * oc.dot(ray.direction);
    float c = oc.dot(oc) - sphereRadius * sphereRadius;

    float discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return std::nullopt;  // No intersection
    }

    // Return closest intersection point
    float t = (-b - std::sqrt(discriminant)) / (2.0f * a);

    if (t < 0) {
        t = (-b + std::sqrt(discriminant)) / (2.0f * a);
    }

    if (t < 0) {
        return std::nullopt;  // Behind camera
    }

    return t;
}

std::optional<RayHit> RayCaster::castRay(
    const Ray& ray,
    const std::vector<float>& atomSpheres,
    const Mat4& modelMatrix
) {
    float closestDistance = std::numeric_limits<float>::max();
    std::optional<RayHit> closestHit;

    // Each atom sphere has 7 floats: x, y, z, radius, r, g, b
    size_t numAtoms = atomSpheres.size() / 7;

    for (size_t i = 0; i < numAtoms; ++i) {
        size_t offset = i * 7;

        // Get atom position (in model space, but since you don't transform, it's world space)
        Vec3 atomCenter(
            atomSpheres[offset + 0],
            atomSpheres[offset + 1],
            atomSpheres[offset + 2]
        );
        float atomRadius = atomSpheres[offset + 3];

        // Test intersection
        auto t = raySphereIntersection(ray, atomCenter, atomRadius);

        if (t && *t < closestDistance) {
            closestDistance = *t;
            Vec3 hitPoint = ray.origin + ray.direction * (*t);
            closestHit = RayHit(i, *t, hitPoint);
        }
    }

    return closestHit;
}