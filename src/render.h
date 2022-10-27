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

// Change the color of the ray according to Phong Shading model
extern bool drawDebugShading;

// Change the maximum allowed ray depth
extern int max_ray_depth;

// The level in the recursion tree of which to show the intersected but unvisited nodes of the BVH
extern bool showUnvisited;
extern int traversalDebugDepth;

extern bool drawDebugSupersamplingRays;
extern int samplesPerPixel;

extern float focalLength;
extern float aperture;
extern int DOFsamples;

// Main rendering function.
void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);

float getRand(float x = 0.0f, float y = 1.0f - std::numeric_limits<float>::epsilon());

glm::vec3 pixelColorDOF(const Scene& scene, const BvhInterface& bvh, Ray& ray, const Features& features, const int rayDepth);

// Get the color of a ray.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth = 1);