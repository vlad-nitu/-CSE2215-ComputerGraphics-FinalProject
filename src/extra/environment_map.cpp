#include "environment_map.h"
#include <render.h>
#include <texture.h>
#include <framework/image.h>
#include <numbers>

Image environment = Image(".\\..\\..\\..\\data\\environment.png");

glm::vec3 getEnvironmentColor(const glm::vec3& lightDirection, const Features& features)
{
    glm::vec3 dir = glm::normalize(lightDirection);

    float u = glm::atan(dir.x, dir.z) / (2 * std::numbers::pi) + 0.5f;
    float v = dir.y * 0.5f + 0.5f;

    glm::vec2 texCoord = glm::vec2 { u, v };

    return acquireTexel(environment, texCoord, features);
}
