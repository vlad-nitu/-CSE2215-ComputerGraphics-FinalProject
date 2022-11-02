#include "environment_map.h"
#include <render.h>
#include <texture.h>
#include <framework/image.h>
#include <numbers>

Image environment = Image(".\\..\\..\\..\\data\\environment.png");
std::vector<Image> cubeMap = {
    // Positive and negative x
    Image(".\\..\\..\\..\\data\\right.png"),
    Image(".\\..\\..\\..\\data\\left.png"),

    // Positive and negative y
    Image(".\\..\\..\\..\\data\\top.png"),
    Image(".\\..\\..\\..\\data\\bottom.png"),

    // Positive and negative z
    Image(".\\..\\..\\..\\data\\front.png"),
    Image(".\\..\\..\\..\\data\\back.png"),
};

glm::vec3 getEnvironmentColor(const glm::vec3& lightDirection, const Features& features)
{
    glm::vec3 dir = glm::normalize(lightDirection);

    float u = glm::atan(dir.x, dir.z) / (2 * std::numbers::pi) + 0.5f;
    float v = dir.y * 0.5f + 0.5f;

    glm::vec2 texCoord = glm::vec2 { u, v };

    return acquireTexel(environment, texCoord, features);
}

int getImageAndCoord(const glm::vec3& lightDirection, float& u, float& v)
{
    float absX = fabs(lightDirection.x);
    float absY = fabs(lightDirection.y);
    float absZ = fabs(lightDirection.z);

    int index;

    int isXPositive = lightDirection.x > 0 ? true : false;
    int isYPositive = lightDirection.y > 0 ? true : false;
    int isZPositive = lightDirection.z > 0 ? true : false;

    float maxAxis, uc, vc;

    // POSITIVE X
    if (isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from +z to -z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = -lightDirection.z;
        vc = lightDirection.y;
        index = 0;
    }
    // NEGATIVE X
    if (!isXPositive && absX >= absY && absX >= absZ) {
        maxAxis = absX;
        uc = lightDirection.z;
        vc = lightDirection.y;
        index = 1;
    }
    // POSITIVE Y
    if (isYPositive && absY >= absX && absY >= absZ) {
        maxAxis = absY;
        uc = lightDirection.x;
        vc = -lightDirection.z;
        index = 2;
    }
    // NEGATIVE Y
    if (!isYPositive && absY >= absX && absY >= absZ) {
        maxAxis = absY;
        uc = lightDirection.x;
        vc = lightDirection.z;
        index = 3;
    }
    // POSITIVE Z
    if (isZPositive && absZ >= absX && absZ >= absY) {
        maxAxis = absZ;
        uc = lightDirection.x;
        vc = lightDirection.y;
        index = 4;
    }
    // NEGATIVE Z
    if (!isZPositive && absZ >= absX && absZ >= absY) {
        maxAxis = absZ;
        uc = -lightDirection.x;
        vc = lightDirection.y;
        index = 5;
    }

    // Convert range from -1 to 1 to 0 to 1
    u = 0.5f * (uc / maxAxis + 1.0f);
    v = 0.5f * (vc / maxAxis + 1.0f);

    return index;
}

glm::vec3 getCubeMapColor(const glm::vec3& lightDirection, const Features& features)
{
    float u;
    float v;
    int index = getImageAndCoord(lightDirection, u, v);

    return acquireTexel(cubeMap[index], glm::vec2 { u, v }, features);
}
