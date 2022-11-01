#include <render.h>
#include "screen.h"
#include <framework/trackball.h>
#include "supersampling.h"

Ray generateSample(const int x, const int y, const int p, const int q, const glm::ivec2& windowResolution, const Trackball& camera)
{
    glm::vec2 samplePosition {
        (float(x + (float(p) + getRand()) / samplesPerPixel) / float(windowResolution.x)) * 2.0f - 1.0f,
        (float(y + (float(q) + getRand()) / samplesPerPixel) / float(windowResolution.y)) * 2.0f - 1.0f
    };

    return camera.generateRay(samplePosition);
}

glm::vec3 pixelResult(const int sampleCount, const int x, const int y, const glm::ivec2& windowResolution, const Trackball& camera, const Scene& scene, const BvhInterface& bvh, const Features& features)
{
    glm::vec3 pixelColor = glm::vec3 { 0 };

    for (int p = 0; p < sampleCount; p++) {
        for (int q = 0; q < sampleCount; q++) {

            pixelColor += getFinalColor(scene, bvh, generateSample(x, y, p, q, windowResolution, camera), features);
        }
    }

    pixelColor /= (sampleCount * sampleCount);
    return pixelColor;
}

void supersampling(const int sampleCount, const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            screen.setPixel(x, y, pixelResult(sampleCount, x, y, windowResolution, camera, scene, bvh, features));
        }
    }
}
