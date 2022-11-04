#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>

#include <glm/gtc/random.hpp> // Ask TA if this is allowed
#include <random> // Ask Ta if this is allowed

// Import in order to perform second visual debug for BVH traversal
#include <bounding_volume_hierarchy.h> // Ask TAs if allowed

#include <extra/environment_map.h>
#include <extra/supersampling.h>

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

bool useSphereEnvironment = false;

int samplesPerPixel = 2; // Sample size per pixel

// Set focal length for depth of field
float focalLength = 4.0f;
float aperture = 0.05f;
int DOFsamples = 1;

// Constant seed used when debugging
// N.B.: Not using a dynamically-computed seed (the random device and the seed are computed only once every run)
// Otherwise, once the project is run, the perturbed samples would be continuously recalculated, and debugging would be difficult.
std::random_device rdGlossyConstant;
unsigned seedGlossyConstant = rdGlossyConstant();

bool glossyConstantSeed = false;
float degreeOfBlur = 0.25f;
int numPerturbedSamples = 50;

bool bloomDebug = false;
float threshold = 0.9f;
int filterSize = 40;
float scalingFactor = 1.0f;

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

                        // Dynamically-computed seed used when rendering
                        std::random_device rdDynamic;
                        unsigned seedDynamic = rdDynamic();

                        // Distinguish between debugging (constant seed) or rendering (dynamically-computed seed)
                        unsigned seed = 0;
                        if (glossyConstantSeed) {
                            seed = seedGlossyConstant;
                        } else {
                            seed = seedDynamic;
                        }

                        // Used for generating uniform random numbers in the interval [0; 1)
                        std::default_random_engine engine(seed);
                        std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

                        // Compute degree of blur -> determined via the corresponding slider
                        // -> the default value is set to 0.25 (= 1 / hitInfo.material.shininess)

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

        if (features.extra.enableEnvironmentMapping) {

            if (useSphereEnvironment)
                Lo += getEnvironmentColor(ray.direction, features);
            else
                Lo += getCubeMapColor(ray.direction, features);

            if (drawDebugShading)
                drawRay(ray, Lo);
            else
                drawRay(ray, glm::vec3 { 1 });

            if (rayDepth == 1) {
                // Implement DOF
                if (features.extra.enableDepthOfField) {
                    Lo /= DOFsamples; // Average the results without main ray
                }
            }

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
        glm::vec3 randomPos = glm::vec3 { getRand(-aperture, aperture), getRand(-aperture, aperture), getRand(-aperture, aperture) };

        glm::vec3 lense = randomPos - (glm::dot(randomPos, ray.direction) / glm::dot(ray.direction, ray.direction)) * ray.direction;
        lense = glm::normalize(lense);

        glm::vec3 randomLense = ray.origin + lense * getRand(-aperture, aperture);

        Ray newRay = { randomLense, glm::normalize(focalPoint - randomLense), std::numeric_limits<float>::max() };

        color += getFinalColor(scene, bvh, newRay, features, rayDepth + 1);
    }

    return color;
}

// Convert the given channel (R, G, or B) to a linear value (between 0 and 1)
// Needed for the relative luminosity formula -> used in bloomFilter()
// Reference used: https://en.wikipedia.org/wiki/SRGB
float toLinearRGB(float channelValue)
{
    if (channelValue > 0.04045f) {

        return powf(((channelValue + 0.055f) / 1.055f), 2.4f);

    } else {

        return channelValue / 12.92f;

    }
}

// Apply box filter (as shown in the slides for lecture 2 of the course)
glm::vec3 boxFilter(std::vector<std::vector<glm::vec3>>& pixelGrid, int x, int y, int givenFilterSize)
{
    glm::vec3 result = glm::vec3(0.0f);

    for (int i = -givenFilterSize; i < (givenFilterSize + 1); i++) {
        for (int j = -givenFilterSize; j < (givenFilterSize + 1); j++) {

            result += pixelGrid[x + i][y + j];

        }
    }

    int denominator = (2 * givenFilterSize + 1) * (2 * givenFilterSize + 1);
    return (result / static_cast<float>(denominator));
}

// Apply bloom filter (only keeping colors with a relative luminance greater than the chosen threshold)
void bloomFilter(glm::ivec2& resolution, std::vector<std::vector<glm::vec3>>& pixelGrid, float givenThreshold, int givenFilterSize)
{
    for (int y = 0; y < resolution.y; y++) {
        for (int x = 0; x < resolution.x; x++) {

            glm::vec3 originalColor = pixelGrid[x + givenFilterSize][y + givenFilterSize];

            // Convert the current RGB values to linear values (between 0 and 1)
            glm::vec3 convertedToLinear = glm::vec3(toLinearRGB(originalColor.x), toLinearRGB(originalColor.y), toLinearRGB(originalColor.z));

            // Calculating the relative luminance of the current RGB values
            // Reference used: https://en.wikipedia.org/wiki/Relative_luminance
            float valueToCheck = glm::dot(convertedToLinear, glm::vec3(0.2126f, 0.7152f, 0.0722f));

            if (valueToCheck > givenThreshold) {

                // If the relative luminance meets the threshold, keep the original color
                continue;

            } else {

                // Otherwise, replace the original color with black [RGB(0, 0, 0)]
                pixelGrid[x + givenFilterSize][y + givenFilterSize] = glm::vec3(0.0f);

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
        std::vector<glm::vec3> originalColors = screen.pixels();

        // Create an extended grid containing the screen's original pixels and a border of length "filterSize"
        // -> The border will come into place once the original colors are stored within the grid
        // -> It surrounds the internal-screen positions from each side (top, right, bottom, left)
        // -> Its purpose is to ensure that box filtering will be possible everywhere across the internal-screen positions
        // -> (particularly those at the edges of the screen's original pixels)
        // -> (there, the box filter would otherwise attempt to access out-of-bounds positions)
        // 
        // Initialize every position (both borders and internal-screen positions) with black [RGB(0, 0, 0)]
        std::vector<std::vector<glm::vec3>> extendedPixelGrid(size_t(windowResolution.x + 2 * filterSize));
        for (int i = 0; i < windowResolution.x + 2 * filterSize; i++) {

            std::vector<glm::vec3> rowWithBorder(size_t(windowResolution.y + 2 * filterSize), glm::vec3(0.0f));
            extendedPixelGrid[i] = rowWithBorder;

        }

        // Store the screen's original pixels within the grid (this way also creating the border)
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {

                // Introduce offsets of "filterSize" to keep the border intact
                extendedPixelGrid[x + filterSize][y + filterSize] = originalColors[screen.indexAt(x, y)];

            }
        }

        // Perform bloom filter (keep only those values with a relative luminosity greater than the chosen threshold)
        bloomFilter(windowResolution, extendedPixelGrid, threshold, filterSize);

        for (int y = filterSize; y < windowResolution.y + filterSize; y++) {
            for (int x = filterSize; x < windowResolution.x + filterSize; x++) {

                if (bloomDebug) {

                    // Visual debug -> display ONLY the bloom filter pixels (nothing else)
                    // N.B.: all other previously-computed information is discarded (this checkbox should ONLY be used for debugging)
                    // -> when the aim is to render a complete image, uncheck the "bloomDebug" box
                    // -> with this separation, the effects of the chosen scalingFactor and filterSize become more easily-interpretable
                    screen.setPixel(x - filterSize, y - filterSize, scalingFactor * boxFilter(extendedPixelGrid, x, y, filterSize));

                } else {

                    // Add the resulting color to the corresponding original pixel (multiplying by the scalingFactor to regulate the brightness)
                    screen.setPixel(x - filterSize, y - filterSize, originalColors[screen.indexAt(x - filterSize, y - filterSize)] + scalingFactor * boxFilter(extendedPixelGrid, x, y, filterSize));

                }

            }
        }
    }
}