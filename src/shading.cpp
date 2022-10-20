#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeDiffuse(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 postion = ray.origin + ray.t * ray.direction;

    glm::vec3 surfaceLightVector = glm::normalize(lightPosition - postion);
    float dotProduct = std::max(0.0f, glm::dot(surfaceLightVector, hitInfo.normal));

    if (features.enableTextureMapping && hitInfo.material.kdTexture) {
        glm::vec3 kd = acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);

        return lightColor * kd * dotProduct;
    } else {
        return lightColor * hitInfo.material.kd * dotProduct;
    }
}

const glm::vec3 computeSpecular(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 position = ray.origin + ray.t * ray.direction;

    if (glm::dot(hitInfo.normal, lightPosition - position) <= 0) {
        return glm::vec3 { 0 };
    } else {
        glm::vec3 lightRay = glm::normalize(position - lightPosition); // From light to point

        glm::vec3 reflection = glm::normalize(lightRay - 2 * glm::dot(hitInfo.normal, lightRay) * hitInfo.normal); // from point to infinity
        glm::vec3 viewDirection = glm::normalize(ray.origin - position); // from point to camera

        float incomingLightCoeff = glm::dot(reflection, viewDirection);
        if (incomingLightCoeff > 0) {
            float intensity = pow(incomingLightCoeff, hitInfo.material.shininess);

            glm::vec3 result = lightColor * hitInfo.material.ks * intensity;

            return result;
        } else {
            return glm::vec3 { 0 };
        }
    }
}

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        glm::vec3 diffuse = computeDiffuse(lightPosition, lightColor, features, ray, hitInfo);
        glm::vec3 specular = computeSpecular(lightPosition, lightColor, features, ray, hitInfo);

        glm::vec3 shading = diffuse + specular;
        
        return shading;
    } else {
        return hitInfo.material.kd;
    }
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    // TODO: implement the reflection ray computation.
    return reflectionRay;
}