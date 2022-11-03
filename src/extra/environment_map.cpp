#include "environment_map.h"
#include <render.h>
#include <texture.h>
#include <framework/image.h>
#include <numbers>
#include <draw.h>

Image environment = Image(".\\..\\..\\..\\data\\environment.png");

bool drawEdgeRays = false;

std::vector<glm::vec3> debugRays = {
    glm::vec3 { 1, 1, 1 }, // 0
    glm::vec3 { 1, 1, -1 }, // 1
    glm::vec3 { 1, -1, 1 }, // 2
    glm::vec3 { 1, -1, -1 }, // 3
    glm::vec3 { -1, 1, 1 }, // 4
    glm::vec3 { -1, 1, -1 }, // 5
    glm::vec3 { -1, -1, 1 }, // 6
    glm::vec3 { -1, -1, -1 }, // 7
};
std::vector<glm::vec4> debugRaysPositions = {
    glm::vec4 {0, 1, 2, 3},
    glm::vec4 {4, 6, 5, 7},
    glm::vec4 {0, 4, 5, 1},
    glm::vec4 {2, 6, 7, 3},
    glm::vec4 {2, 6, 4, 0},
    glm::vec4 {1, 5, 7, 3},
};
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
    // Compute in order to compare the power of each axis
    float absX = fabs(lightDirection.x); 
    float absY = fabs(lightDirection.y);
    float absZ = fabs(lightDirection.z);

    int index;

    int isXPositive = lightDirection.x > 0 ? true : false;
    int isYPositive = lightDirection.y > 0 ? true : false;
    int isZPositive = lightDirection.z > 0 ? true : false;

    float maxAxis, uc, vc;

    // POSITIVE X and it it the most powerful in absolute term
    if (isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from +z to -z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = -lightDirection.z;
        vc = lightDirection.y;
        index = 0;
    }
    // NEGATIVE X and it it the most powerful in absolute term
    if (!isXPositive && absX >= absY && absX >= absZ) {
        maxAxis = absX;
        uc = lightDirection.z;
        vc = lightDirection.y;
        index = 1;
    }
    // POSITIVE Y and it it the most powerful in absolute term
    if (isYPositive && absY >= absX && absY >= absZ) {
        maxAxis = absY;
        uc = lightDirection.x;
        vc = -lightDirection.z;
        index = 2;
    }
    // NEGATIVE Y and it it the most powerful in absolute term
    if (!isYPositive && absY >= absX && absY >= absZ) {
        maxAxis = absY;
        uc = lightDirection.x;
        vc = lightDirection.z;
        index = 3;
    }
    // POSITIVE Z and it it the most powerful in absolute term
    if (isZPositive && absZ >= absX && absZ >= absY) {
        maxAxis = absZ;
        uc = lightDirection.x;
        vc = lightDirection.y;
        index = 4;
    }
    // NEGATIVE Z and it it the most powerful in absolute term
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

    if (drawEdgeRays) {
        drawRay({ glm::vec3 { 0 }, debugRays[debugRaysPositions[index].x] }, glm::vec3 { 0, 1, 0.76 });
        drawRay({ glm::vec3 { 0 }, debugRays[debugRaysPositions[index].y] }, glm::vec3 { 0, 1, 0.76 });
        drawRay({ glm::vec3 { 0 }, debugRays[debugRaysPositions[index].z] }, glm::vec3 { 0, 1, 0.76 });
        drawRay({ glm::vec3 { 0 }, debugRays[debugRaysPositions[index].w] }, glm::vec3 { 0, 1, 0.76 });

        drawRay({ glm::vec3 { 0 }, lightDirection, 1.0f }, glm::vec3 { 1, 1, 0 });
    }

    return acquireTexel(cubeMap[index], glm::vec2 { u, v }, features);
}
