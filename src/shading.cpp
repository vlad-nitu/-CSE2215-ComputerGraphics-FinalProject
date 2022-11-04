#include "draw.h"
#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

bool drawReflectionDebug = false;

std::vector<Image> createImages(const Image& image) { 

    std::vector<Image> images{};
    images.push_back(image); 

    for (int lvl = 1; lvl <= mipmap_max_depth; ++lvl) { 

        const Image& prev_image = images.back(); 

        int prev_w = prev_image.width;
        int prev_h = prev_image.height;

        int new_w = prev_w / 2;
        int new_h = prev_h / 2;


        std::vector<glm::vec3> new_pixels {}; 

        int x, y;
        x = y = 0;
        int grid_size = prev_h * prev_w;

        while (y < grid_size) {

            glm::vec3 interpolated_pixel = (prev_image.pixels[x + y] + prev_image.pixels[x + y + 1] + prev_image.pixels[x + y + prev_w] + prev_image.pixels[x + y + prev_w + 1]) / 4.0f;
            new_pixels.push_back(interpolated_pixel);

            if (x < prev_w)
            {
                x += 2;
                if (x == prev_w)
                {
                    x = 0; y += 2 * prev_w;
                }
            }
            else{
                x = 0; y += 2 * prev_w;
            }
        }

        Image new_image = Image(new_w, new_h, new_pixels);
        images.push_back(new_image);
    }
    return images;
}

const glm::vec3 computeDiffuse(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    
    glm::vec3 postion = ray.origin + ray.t * ray.direction;

    glm::vec3 surfaceLightVector = glm::normalize(lightPosition - postion);
    float dotProduct = std::max(0.0f, glm::dot(surfaceLightVector, hitInfo.normal));

    if (features.enableTextureMapping && features.extra.enableMipmapTextureFiltering && hitInfo.material.kdTexture){

            Image& img = *hitInfo.material.kdTexture;
             int level;
             if (ray.t <= 1.0f) // avoid log2 from 0
                 level = 0;
             else
                 level = log2(ray.t);

            mipmap_max_depth = std::log2(img.height);

            if (map.find(img) == map.end())
                map[img] = createImages(img);

            if (level > mipmap_max_depth)
                level = mipmap_max_depth;

            float alpha;
            if (level == mipmap_max_depth)
            {
                // only use the last mipmap level, the one that contains 1x1 pixels
                level --;
                alpha = 1; // only use the last mipmap level, the one that contains 1x1 pixels
            }
            else
                alpha = (ray.t - pow(2, level)) / (pow(2, level + 1) - pow(2, level));

            if (features.enableTextureMapping && hitInfo.material.kdTexture)
            {
                glm::vec3 kd;
                
                if (features.extra.enableBilinearTextureFiltering)
                {
                    kd = (1 - alpha) * bilinearInterpolation(map[img][level], hitInfo.texCoord, features, level, ray) + alpha * bilinearInterpolation(map[img][level + 1], hitInfo.texCoord, features, level, ray) ;
                }
                else
                {
                    kd = (1 - alpha) * acquireTexel(map[img][level], hitInfo.texCoord, features, level, ray) + alpha * acquireTexel(map[img][level + 1], hitInfo.texCoord, features, level, ray) ;
                }

                return lightColor * kd * dotProduct;
            }
            else 
                return lightColor * hitInfo.material.kd * dotProduct;
    }

    if (features.enableTextureMapping && hitInfo.material.kdTexture) {
        glm::vec3 kd; 
        if (features.extra.enableBilinearTextureFiltering)
            kd = bilinearInterpolation(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
        else
            kd = acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);

        return lightColor * kd * dotProduct;
    } else {
        return lightColor * hitInfo.material.kd * dotProduct;
    }
}

const glm::vec3 computeSpecular(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray,
    HitInfo hitInfo)
{
    glm::vec3 position = ray.origin + ray.t * ray.direction;

    if (glm::dot(hitInfo.normal, lightPosition - position) <= 0) {
        return glm::vec3 { 0 };
    } else {
        glm::vec3 lightRay = glm::normalize(position - lightPosition); // From light to point

        glm::vec3 reflection = glm::normalize(
            lightRay - 2 * glm::dot(hitInfo.normal, lightRay) * hitInfo.normal); // from point to infinity
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

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray,
    HitInfo hitInfo)
{
    if (features.enableShading) {
        glm::vec3 diffuse = computeDiffuse(lightPosition, lightColor, features, ray, hitInfo);
        glm::vec3 specular = computeSpecular(lightPosition, lightColor, features, ray, hitInfo);

        glm::vec3 shading = diffuse + specular;

        return shading;
    } else {
        if (features.enableTextureMapping && hitInfo.material.kdTexture) {

             if (features.extra.enableMipmapTextureFiltering) {

             Image& img = *hitInfo.material.kdTexture;
             int level;
             if (ray.t <= 1.0f) // avoid log2 from 0
                 level = 0;
             else
                 level = log2(ray.t);

             mipmap_max_depth = std::log2(img.height);

            if (map.find(img) == map.end())
                map[img] = createImages(img);

             if (level > mipmap_max_depth)
                level = mipmap_max_depth;

             float alpha;
             if (level == mipmap_max_depth)
             {
                 // only use the last mipmap level, the one that contains 1x1 pixels
                 level --;
                 alpha = 1;
             }
             else
                 alpha = (ray.t - pow(2, level)) / (pow(2, level + 1) - pow(2, level));

            if (features.extra.enableBilinearTextureFiltering)
                return (1 - alpha) * bilinearInterpolation(map[img][level], hitInfo.texCoord, features, level, ray) + alpha * bilinearInterpolation(map[img][level + 1], hitInfo.texCoord, features, level, ray) ;
            else
                return (1 - alpha) * acquireTexel(map[img][level], hitInfo.texCoord, features, level, ray) + alpha * acquireTexel(map[img][level + 1], hitInfo.texCoord, features, level, ray) ;
             }

            if (features.extra.enableBilinearTextureFiltering)
                return bilinearInterpolation(*hitInfo.material.kdTexture, hitInfo.texCoord, features); 
            else
                return acquireTexel(*hitInfo.material.kdTexture, hitInfo.texCoord, features);
        } else {
            return hitInfo.material.kd;
        }
    }
}

const Ray computeReflectionRay(Ray ray, HitInfo hitInfo)
{
    glm::vec3 point = ray.origin + ray.t * ray.direction;

    glm::vec3 N = hitInfo.normal;
    glm::vec3 L = glm::normalize(ray.direction);

    // Reflect L over N -> R = L - 2 (L*N) * N
    glm::vec3 R = L - 2.0f * glm::dot(L, N) * N;
    glm::vec3 normalized_R = glm::normalize(R);

    const float ERR = 1e-5;
    Ray reflectionRay { point + ERR * normalized_R, normalized_R, std::numeric_limits<float>::max() };
    if (drawReflectionDebug) {
        drawRay({ reflectionRay.origin, reflectionRay.direction, 0.3f }, glm::vec3 { 0, 0, 1 });
    }
    return reflectionRay;
}

const Ray computeRefractedRay(const Ray& ray)
{
    const float OFFSET = 1e-5;
    glm::vec3 p = ray.origin + ray.t * ray.direction;

    return {p + OFFSET * ray.direction, ray.direction, std::numeric_limits<float>::max()};
}

