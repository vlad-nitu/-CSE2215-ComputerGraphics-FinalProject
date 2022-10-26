#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <random>

bool drawShadowRayDebug = false;
// Constant seed used for uniform random numbers (for sampling the line-segment and parallelogram light sources)
// N.B.: Not using a random seed (the random device and the seed are computed only once every run)
// Otherwise, once the project is run, the shadow rays would be continuously recalculated, and no static shadows would be produced.
std::random_device rd;
unsigned seed = rd();

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    // Extracting the uniform random number (alpha) -> changing method signatures is not allowed
    float alpha = position.x;

    // Computing the position coordinates (somewhere along the line-segment)
    float coordX = alpha * segmentLight.endpoint1.x + (1.0f - alpha) * segmentLight.endpoint0.x;
    float coordY = alpha * segmentLight.endpoint1.y + (1.0f - alpha) * segmentLight.endpoint0.y;
    float coordZ = alpha * segmentLight.endpoint1.z + (1.0f - alpha) * segmentLight.endpoint0.z;
    position = glm::vec3(coordX, coordY, coordZ);

    // Computing the color (via linear interpolation)
    color = (1.0f - alpha) * segmentLight.color0 + alpha * segmentLight.color1;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
// Reference used: Fundamentals of Computer Graphics (Fourth Edition), Chapter 13, Section 13.4.2, pp. 331-332
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    // Extracting the uniform random numbers (alpha, beta) -> changing method signatures is not allowed
    float alpha = position.x;
    float beta = color.x;

    // Computing the position coordinates (somewhere within the parallelogram)
    position = parallelogramLight.v0 + alpha * parallelogramLight.edge01 + beta * parallelogramLight.edge02;

    // Computing the color (via bilinear interpolation)
    color = alpha * beta * parallelogramLight.color3 + (1.0f - alpha) * beta * parallelogramLight.color2
        + (1.0f - alpha) * (1.0f - beta) * parallelogramLight.color0 + alpha * (1.0f - beta) * parallelogramLight.color1;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
// Reference used: Fundamentals of Computer Graphics (Fourth Edition), Chapter 4, Section 4.7, pp. 86-87
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 currentPoint = ray.origin + ray.t * ray.direction;
    glm::vec3 directionPointToLight = glm::normalize(samplePos - currentPoint);
    float tLight = glm::distance(samplePos, currentPoint);

    Ray pointTowardsLight { currentPoint + 0.0001f * directionPointToLight, directionPointToLight, tLight };
    if (bvh.intersect(pointTowardsLight, hitInfo, features)) {
        if (pointTowardsLight.t > tLight || fabs(pointTowardsLight.t - tLight) < 0.0001f) {
            if (drawShadowRayDebug)
                drawRay(pointTowardsLight, debugColor);
            return 1.0f;
        } else {
            if (drawShadowRayDebug)
                drawRay(pointTowardsLight, glm::vec3(1, 0, 0));
            return 0.0f;
        }
    } else {
        if (drawShadowRayDebug)
            drawRay(pointTowardsLight, debugColor);
        return 1.0f;
    }
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.

        glm::vec3 result = glm::vec3 { 0.0f };

        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);

                if (features.enableHardShadow) {
                    float lightContribution = testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo);

                    result += lightContribution * computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                } else
                    result += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);

            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);

                if (features.enableSoftShadow) {

                    const int SAMPLE_COUNT = 20;

                    glm::vec3 position = glm::vec3 { 0 };
                    glm::vec3 color = glm::vec3 { 0 };

                    // Generating a uniform random number in the interval [0; 1)
                    std::default_random_engine engine(seed);
                    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

                    for (int i = 0; i < SAMPLE_COUNT; i++) {
                        glm::vec3 samplePosition = glm::vec3 { uniform(engine) };
                        glm::vec3 sampleColor = glm::vec3 { 0 };

                        sampleSegmentLight(segmentLight, samplePosition, sampleColor);

                        // Accumulating positions
                        position += samplePosition;

                        // Accumulating colors -> black [RGB(0, 0, 0)] in case the sample position is not visible
                        float isVisible = testVisibilityLightSample(samplePosition, sampleColor, bvh, features, ray, hitInfo);
                        color += isVisible * sampleColor;
                    }
                    // Computing averages for position and color
                    position /= SAMPLE_COUNT;
                    color /= SAMPLE_COUNT;

                    result += computeShading(position, color, features, ray, hitInfo);
                }

            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);

                if (features.enableSoftShadow) {

                    const int SAMPLE_COUNT = 40;

                    glm::vec3 position = glm::vec3 { 0 };
                    glm::vec3 color = glm::vec3 { 0 };

                    // Generating two uniform random numbers in the interval [0; 1)
                    std::default_random_engine engine(seed);
                    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

                    for (int i = 0; i < SAMPLE_COUNT; i++) {
                        glm::vec3 samplePosition = glm::vec3 { uniform(engine) };
                        glm::vec3 sampleColor = glm::vec3 { uniform(engine) };

                        sampleParallelogramLight(parallelogramLight, samplePosition, sampleColor);

                        // Accumulating positions
                        position += samplePosition;

                        // Accumulating colors -> black [RGB(0, 0, 0)] in case the sample position is not visible
                        float isVisible = testVisibilityLightSample(samplePosition, sampleColor, bvh, features, ray, hitInfo);
                        color += isVisible * sampleColor;
                    }
                    // Computing averages for position and color
                    position /= SAMPLE_COUNT;
                    color /= SAMPLE_COUNT;

                    result += computeShading(position, color, features, ray, hitInfo);
                }
            }
        }
        // return hitInfo.material.kd;
        return result;

    } else {
        return computeShading(glm::vec3 { 0 }, glm::vec3 { 0 }, features, ray, hitInfo);
    }
}
