#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>

#include <glm/gtc/random.hpp> // Ask TA if this is allowed
#include <random> // Ask Ta if this is allowed

#ifdef NDEBUG
#include <omp.h>
#endif

bool drawDebugShading = false;
int ray_depth = 0;

bool drawDebugSupersamplingRays = false;

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);


        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            // TODO: put your own implementation of recursive ray tracing here.
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

/// <summary>
/// Generates a random float value between x and y inclusive. By default returns in [0, 1)
/// </summary>
/// <param name="x"> Lower bound, by default 0 </param>
/// <param name="y"> Upper bound, by default 1 - max_epsilon </param>
/// <returns> Returns a random float </returns>
float getRand(float x, float y)
{
    //return glm::linearRand(x, y);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(x, y);

    return dis(gen);
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
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

            // Check if we need to turn on multiple rays per pixel
            if (features.extra.enableMultipleRaysPerPixel) {
                const int n = 3; // The number of total samples is represented by n * n

                glm::vec3 pixelColor = glm::vec3 { 0 };

                /*
                * Implementaion taken from:
                * 
                * Fundamentals of Computer Graphics, 4th Edition, Steve Marschner and Peter Shirley
                * Chapter 13.4.1
                * 
                * In order to perform irregular sampling we need to introduce some sort of randomness into our computations
                * However full randomness can create some problems, such as random patterns. This is why we have chosen to implement
                * an in-between algorithm. We subdive the pixle into n^2 smaller pixels and for each one of them we cast a random ray.
                * 
                * This method retains the random property while making sure that no clusters or patterns arrise.
                */
                for (int p = 0; p < n; p++) {
                    for (int q = 0; q < n; q++) {

                        // Compute the position inside the pixel where this ray should go through
                        glm::vec2 samplePosition{
                            (float(x + (float(p) + getRand()) / n) / float(windowResolution.x)) * 2.0f - 1.0f,
                            (float(y + (float(q) + getRand()) / n) / float(windowResolution.y)) * 2.0f - 1.0f
                        };
                        const Ray sampleRay = camera.generateRay(samplePosition); // Compute the ray

                        pixelColor += getFinalColor(scene, bvh, sampleRay, features);
                    }
                }

                // Average the values
                pixelColor /= n * n;
                screen.setPixel(x, y, pixelColor);

            } else {
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
            }
        }
    }
}