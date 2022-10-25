#pragma once

#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include "common.h"
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()

// Forward declarations.
struct Image;

// Given an image and a texture coordinate, return the corresponding texel.
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features, int level);
int getMipMapLevel(int h);
std::vector<Image> createImages(const Image& image);