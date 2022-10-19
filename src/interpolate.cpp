#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 normal_largeTriangle = glm::cross((v1 - v0), (v2 - v0));
    glm::vec3 normal_smallTriangleOppositeV0 = glm::cross((v2 - v1), (p - v1));
    glm::vec3 normal_smallTriangleOppositeV1 = glm::cross((v0 - v2), (p - v2));
    glm::vec3 normal_smallTriangleOppositeV2 = glm::cross((v1 - v0), (p - v0));

    float denominator = glm::dot(normal_largeTriangle, normal_largeTriangle);
    float alpha = glm::dot(normal_largeTriangle, normal_smallTriangleOppositeV0) / denominator;
    float beta = glm::dot(normal_largeTriangle, normal_smallTriangleOppositeV1) / denominator;
    float gamma = glm::dot(normal_largeTriangle, normal_smallTriangleOppositeV2) / denominator;

    return glm::vec3(alpha, beta, gamma);
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    glm::vec3 interpolatedNormal = barycentricCoord.x * n0 + barycentricCoord.y * n1 + barycentricCoord.z * n2;
    return glm::normalize(interpolatedNormal);
}

glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return t0 * barycentricCoord.x + t1 * barycentricCoord.y + t2 * barycentricCoord.z;
}
