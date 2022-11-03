#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>

#include <glm/gtc/random.hpp> // Ask TA if this is allowed
#include <random> // Ask Ta if this is allowed

// Import in order to perform second visual debug for BVH traversal
#include <bounding_volume_hierarchy.h> // Ask TAs if allowed

#include <supersampling.h>

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

int samplesPerPixel = 2; // Sample size per pixel

// Set focal length for depth of field
float focalLength = 4.0f;
float aperture = 0.05f;
int DOFsamples = 1;

int numPerturbedSamples = 50;

bool bloomDebug = false;

// Extracted functionality for normal recursive raytracing (used for computing the color of perturbed samples)
glm::vec3 perturbedSampleColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    glm::vec3 Lo = glm::vec3 { 0 };

    if (bvh.intersect(ray, hitInfo, features)) {
        Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (rayDepth <= max_ray_depth) {
            if (hitInfo.material.ks != glm::vec3 { 0, 0, 0 }) {
                Ray reflection = computeReflectionRay(ray, hitInfo);

                Lo += hitInfo.material.ks * perturbedSampleColor(scene, bvh, reflection, features, rayDepth + 1);
            }
        }
    }

    return Lo;
}

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    glm::vec3 Lo = glm::vec3 { 0 };

    if (showUnvisited) {
        if (rayDepth == traversalDebugDepth)
            drawUnvisited = true;
        else
            drawUnvisited = false;
    } else {
        drawUnvisited = false;
    }

    // Check if the rays come from the camera
    if (rayDepth == 1) {
        // Implement DOF
        if (features.extra.enableDepthOfField) {

            // Also trace the random focal point rays
            Lo += pixelColorDOF(scene, bvh, ray, features, rayDepth);
        }
    }

    if (bvh.intersect(ray, hitInfo, features)) {

        Lo += computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {

            // Check if we reached the maximum depth
            if (rayDepth <= max_ray_depth) {

                // Check if the transparency feature is enabled
                if (features.extra.enableTransparency) {

                    float transparency = hitInfo.material.transparency;

                    // Check if material can reflect light
                    if (hitInfo.material.ks != glm::vec3 { 0, 0, 0 }) {
                        Ray reflection = computeReflectionRay(ray, hitInfo);

                        float angle = glm::dot(glm::normalize(ray.direction), glm::normalize(reflection.direction));
                        Lo += transparency * hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
                    }

                    // Check if material is transparent
                    if (transparency < 1.0f) {
                        Ray refraction = computeRefractedRay(ray);

                        Lo *= transparency;

                        Lo += (1.0f - transparency) * getFinalColor(scene, bvh, refraction, features, rayDepth + 1);
                    }

                } else if (features.extra.enableGlossyReflection) {
                    // Check if glossy reflections are enabled

                    // Check if light can be reflected off the material
                    if (hitInfo.material.ks != glm::vec3 { 0, 0, 0 }) {
                        // Compute primary reflection ray and its direction (before perturbing)
                        Ray reflection = computeReflectionRay(ray, hitInfo);
                        glm::vec3 reflectionDirection = reflection.direction;

                        // Compute a vector that is NOT collinear with the reflection direction
                        // Reference used: Fundamentals of Computer Graphics (Fourth Edition), Chapter 2, Section 2.4.6, p. 28
                        float minimumComponent = std::fmin(fabs(reflectionDirection.x), std::fmin(fabs(reflectionDirection.y), fabs(reflectionDirection.z)));
                        glm::vec3 T { 0, 0, 0 };
                        if (fabs(fabs(reflectionDirection.x) - minimumComponent) < 0.0001f) {
                            T = glm::vec3 { 1, reflectionDirection.y, reflectionDirection.z };
                        } else if (fabs(fabs(reflectionDirection.y) - minimumComponent) < 0.0001f) {
                            T = glm::vec3 { reflectionDirection.x, 1, reflectionDirection.z };
                        } else {
                            T = glm::vec3 { reflectionDirection.x, reflectionDirection.y, 1 };
                        }

                        // Compute alligned orthonormal basis -> obtaining three orthonormal vectors (reflectionDirection, U, V)
                        glm::vec3 U = glm::normalize(glm::cross(T, reflectionDirection));
                        glm::vec3 V = glm::cross(reflectionDirection, U);

                        // Used for generating uniform random numbers in the interval [0; 1)
                        std::random_device rd;
                        std::default_random_engine engine(rd());
                        std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

                        // Compute degree of blur -> currently, the inverse of the material shininess
                        // -> the value is 0.25 (as the shininess of the parallelepiped mirror is 4)
                        // Question -> should this be a slider instead?
                        // -> a range of 0.1 to 0.9 produces good images (of increasing blur)
                        float degreeOfBlur = 1 / hitInfo.material.shininess;

                        // Aggregation vector (later used for averaging the colors from all perturbed samples)
                        glm::vec3 aggregatedColors { 0, 0, 0 };

                        for (int i = 0; i < numPerturbedSamples; i++) {
                            float num1 = uniform(engine);
                            float num2 = uniform(engine);

                            // Compute two coefficients (to determine the direction of the current perturbed sample)
                            // Reference used: Fundamentals of Computer Graphics (Fourth Edition), Chapter 13, Section 13.4.4, pp. 333-334
                            float u = -0.5f * degreeOfBlur + num1 * degreeOfBlur;
                            float v = -0.5f * degreeOfBlur + num2 * degreeOfBlur;

                            // Compute the direction of the current perturbed sample and the ray itself
                            glm::vec3 currentDirection = glm::normalize(reflectionDirection + u * U + v * V);
                            Ray perturbedReflection { reflection.origin, currentDirection, std::numeric_limits<float>::max() };

                            // Compute cosine of the angle between the primary reflection ray and the current perturbed sample
                            float angle = glm::dot(reflectionDirection, currentDirection);

                            // Compute the color of the current perturbed sample (using the recursive helper method)
                            // Take sharpness into consideration (via the shininess-exponentiation)
                            // Add the current color to the accumulator
                            glm::vec3 currentColor = powf(angle, hitInfo.material.shininess) * hitInfo.material.ks * perturbedSampleColor(scene, bvh, perturbedReflection, features, rayDepth);
                            aggregatedColors += currentColor;

                            // Visual debug -> draw the perturbed samples (alongside the primary reflection ray colored in blue)
                            // Note -> continuously recalculated due to the dynamically-computed seed
                            if (drawReflectionDebug) {
                                drawRay({ perturbedReflection.origin, perturbedReflection.direction, 0.3f }, currentColor);
                            }
                        }

                        // Add the average color (from all perturbed samples) to Lo
                        Lo += (aggregatedColors / static_cast<float>(numPerturbedSamples));
                    }

                } else {

                    // Perform normal recursive raytracing
                    if (hitInfo.material.ks != glm::vec3 { 0, 0, 0 }) {
                        Ray reflection = computeReflectionRay(ray, hitInfo);

                        float angle = glm::dot(glm::normalize(ray.direction), glm::normalize(reflection.direction));
                        Lo += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
                    }
                }
            }
        }

        // Check if the rays come from the camera
        if (rayDepth == 1) {
            if (features.extra.enableDepthOfField) {
                Lo /= (DOFsamples + 1); // Average the results
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

        // Check if the rays come from the camera
        if (rayDepth == 1) {
            // Implement DOF
            if (features.extra.enableDepthOfField) {
                Lo /= DOFsamples; // Average the results without main ray
            }
        }

        return Lo;
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
    // return glm::linearRand(x, y);

    // Implementation taken from https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(x, y);

    return dis(gen);
}

glm::vec3 pixelColorDOF(const Scene& scene, const BvhInterface& bvh, Ray& ray, const Features& features, int rayDepth)
{
    glm::vec3 color = glm::vec3 { 0 };
    glm::vec3 focalPoint = ray.origin + focalLength * ray.direction;

    for (int i = 0; i < DOFsamples; i++) {
        glm::vec3 randomLense = ray.origin + glm::vec3 { getRand(-aperture, aperture), getRand(-aperture, aperture), getRand(-aperture, aperture) };

        Ray newRay = { randomLense, glm::normalize(focalPoint - randomLense), std::numeric_limits<float>::max() };

        color += getFinalColor(scene, bvh, newRay, features, rayDepth + 1);
    }

    return color;
}

glm::vec3 boxFilter(std::vector<std::vector<glm::vec3>>& pixelMatrix, int x, int y, int sizeProvided)
{
    int filterSize = std::max(1, sizeProvided);
    glm::vec3 result = glm::vec3(0.0f);

    for (int i = -filterSize; i < (filterSize + 1); i++) {
        for (int j = -filterSize; j < (filterSize + 1); j++) {
            result += pixelMatrix[x + i][y + j];
        }
    }

    int denominator = (2 * filterSize + 1) * (2 * filterSize + 1);
    return (result / static_cast<float>(denominator));
}

void bloomFilter(glm::ivec2& resolution, std::vector<std::vector<glm::vec3>>& pixelMatrix, float threshold) {
    for (int y = 0; y < resolution.y; y++) {
        for (int x = 0; x < resolution.x; x++) {

            glm::vec3 originalColor = pixelMatrix[x][y];
            float valueToCheck = glm::dot(originalColor, glm::vec3(0.2126f, 0.7152f, 0.0722f));

            if (valueToCheck > threshold) {
                continue;
            } else {
                pixelMatrix[x][y] = glm::vec3(0.0f);
            }

        }
    }
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode

    if (features.extra.enableMultipleRaysPerPixel) {
        supersampling(samplesPerPixel, scene, camera, bvh, screen, features);
    } else {

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

    if (features.extra.enableBloomEffect) {
        float threshold = 0.9f;
        int filterSize = 15;
        float intensity = 5.0f;

        std::vector<glm::vec3> originalColors = screen.pixels();

        std::vector<std::vector<glm::vec3>> pixelMatrix(windowResolution.x);
        for (int i = 0; i < windowResolution.x; i++) {
            std::vector<glm::vec3> row(windowResolution.y);
            pixelMatrix.push_back(row);
        }

        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                pixelMatrix[x].push_back(originalColors[screen.indexAt(x, y)]);
            }
        }

        bloomFilter(windowResolution, pixelMatrix, threshold);

        for (int y = filterSize; y < windowResolution.y - filterSize; y++) {
            for (int x = filterSize; x < windowResolution.x - filterSize; x++) {

                if (bloomDebug) {
                    screen.setPixel(x, y, intensity * boxFilter(pixelMatrix, x, y, filterSize));
                } else {
                    screen.setPixel(x, y, originalColors[screen.indexAt(x, y)] + intensity * boxFilter(pixelMatrix, x, y, filterSize));
                }

            }
        }
    }
}