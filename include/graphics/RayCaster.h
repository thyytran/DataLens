#pragma once
#ifndef RAYCASTER_H
#define RAYCASTER_H

#include "math/Vec.h"
#include "math/Mat.h"
#include "graphics/Window.h"
#include "Camera.h"
#include "bio/MoleculeData.h"
#include <optional>
#include <limits>

struct Ray {
    Vec3 origin;
    Vec3 direction;

    Ray(const Vec3& origin, const Vec3& direction)
        : origin(origin), direction(direction) {
    }
};

struct RayHit {
    size_t atomIndex;
    float distance;
    Vec3 hitPoint;

    RayHit(size_t atomIndex, float distance, const Vec3& hitPoint)
        : atomIndex(atomIndex), distance(distance), hitPoint(hitPoint) {
    }
};

class RayCaster {
public:
    // Creates a ray from screen coordinates
    static Ray screenToRay(
        double mouseX, double mouseY,
        const Window* window,
        const Camera* camera
    );

    // Tests ray intersection with a sphere
    static std::optional<float> raySphereIntersection(
        const Ray& ray,
        const Vec3& sphereCenter,
        float sphereRadius
    );

    // Finds the closest atom hit by the ray
    static std::optional<RayHit> castRay(
        const Ray& ray,
        const std::vector<float>& atomSpheres,
        const Mat4& modelMatrix
    );

private:
    // Unprojection helper
    static Vec3 unproject(
        const Vec3& screenPos,
        const Mat4& viewMatrix,
        const Mat4& projectionMatrix,
        float screenWidth,
        float screenHeight
    );
};

#endif