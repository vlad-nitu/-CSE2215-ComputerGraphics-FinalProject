#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>

// Forward declarations.
struct Scene;
class Screen;
class Trackball;
class BvhInterface;
struct Features;

extern bool drawDebugShading;
extern int ray_depth;

extern bool drawDebugSupersamplingRays;

extern int samplesPerPixel;

// Main rendering function.
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);

float getRand(float x = 0.0f, float y = 1.0f - std::numeric_limits<float>::epsilon());

// Get the color of a ray.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth = ray_depth);