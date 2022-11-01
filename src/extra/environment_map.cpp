#include "environment_map.h"
#include <render.h>
#include <texture.h>
#include <framework/image.h>
#include <numbers>

Image environment = Image("");

glm::vec3 getEnvironmentColor(const glm::vec3& lightDirection, const Features& features)
{
    float u = glm::atan(lightDirection.x, lightDirection.z) / (2 * std::numbers::pi) + 0.5f;
    float v = lightDirection.y * 0.5f + 0.5f;

    glm::vec2 texCoord = glm::vec2 { u, v };

    return acquireTexel(environment, texCoord, features);
}
