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

Ray generateSample(const int x, const int y, const int p, const int q, const glm::ivec2& windowResolution, const Trackball& camera);

glm::vec3 pixelResult(const int sampleCount, const int x, const int y, const glm::ivec2& windowResolution, const Trackball& camera, const Scene& scene, const BvhInterface& bvh, const Features& features);

void supersampling(const int sampleCount, const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);