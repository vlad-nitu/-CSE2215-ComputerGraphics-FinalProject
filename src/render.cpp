#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>

#ifdef NDEBUG
#include <omp.h>
#endif

glm::vec3 getFinalColor(const Scene &scene, const BvhInterface &bvh, Ray ray, const Features &features, int rayDepth) {
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);


        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);

            if (hitInfo.material.ks != glm::vec3{0, 0, 0} && rayDepth != 0) {
                float angle = glm::dot(glm::normalize(ray.direction),  glm::normalize(reflection.direction));
                Lo = Lo + hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth - 1) * glm::cos(angle);
            }
            // drawRay(reflection, Lo); // debug only atm
        }
        // Draw a debug ray with the color returned from the shading.
        drawRay(ray, Lo);

        // Set the color of the pixel.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

void renderRayTracing(const Scene &scene, const Trackball &camera, const BvhInterface &bvh, Screen &screen,
                      const Features &features) {
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos{
                    float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                    float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }
}