#pragma once

#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include "framework/ray.h"
#include "common.h"
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <unordered_map>
DISABLE_WARNINGS_POP()


// Forward declarations.
struct Image;

extern bool drawMipMapDebug;
namespace std {
    template<>
    struct hash<Image> {
        inline size_t operator()(const Image& img) const {
            std::size_t w_hash = std::hash<int>()(img.width);
            std::size_t h_hash = std::hash<int>()(img.height);
            std::size_t px_hash = std::hash<size_t>()(img.pixels.size());

            return w_hash ^ h_hash ^ px_hash;
        }
    };
}
extern std::unordered_map<Image, std::vector<Image> > map; 
extern int mipmap_max_depth;

// Given an image and a texture coordinate, return the corresponding texel.
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features);
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features, int level, Ray& ray);