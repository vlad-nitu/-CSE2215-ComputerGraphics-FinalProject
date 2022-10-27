#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>

// Import in order to perform second visual debug for BVH traversal
#include <bounding_volume_hierarchy.h> // Ask TAs if allowed

#ifdef NDEBUG
#include <omp.h>
#endif

// Change the color of the ray according to Phong Shading model
bool drawDebugShading = false;

// Change the maximul allowed ray depth
int max_ray_depth = 1;

// The level in the recursion tree of which to show the intersected but unvisited nodes of the BVH
bool showUnvisited = false;
int traversalDebugDepth = 1;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;

    if (showUnvisited) {
        if (rayDepth == traversalDebugDepth)
            drawUnvisited = true;
        else
            drawUnvisited = false;
    } else {
        drawUnvisited = false;
    }

    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);

            if (hitInfo.material.ks != glm::vec3 { 0, 0, 0 } && rayDepth <= max_ray_depth) {
                float angle = glm::dot(glm::normalize(ray.direction), glm::normalize(reflection.direction));
                Lo = Lo + hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
            }
        }

        // Draw a debug ray with the color returned from the shading.
        if (drawDebugShading) {
            drawRay(ray, Lo);
        } else {
            drawRay(ray, glm::vec3 { 1 });
        }

        // Set the color of the pixel.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen,
    const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }
}