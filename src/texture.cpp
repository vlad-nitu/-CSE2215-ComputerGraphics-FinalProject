#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    float u = texCoord.x;
    float v = texCoord.y;

    /*
     * We clamp values to make sure we don't use the wrong pixels for edges or outside normal coordinates
     *
     * This moves the values in the interval [0,1) so that they can be assigned a pixel
     */
    u = std::min(1.0f - 1e-6f, u);
    u = std::max(0.0f, u);

    v = std::min(1.0f - 1e-6f, v);
    v = std::max(0.0f, v);

    int col = u * image.width; // Convert to int so values are rounded down representing the line
    int row = v * image.height; // Convert to int so values are rounded down representing the column

    /*
    * Images are stored top to bottom but coordinates are bottom to top
    * 
    * We need to adjust for this
    * 
    * One is substracted since image.height has values [1, x] and row [0, x-1]
    */
    row = (image.height - row - 1);

    return image.pixels[row * image.width + col];
}

void clampTexel(glm::vec2& center, const Image& image) {
    if (center.x == image.width)
        center.x --;
    if (center.y == image.height)
        center.y --; 
}

glm::vec3 bilinearInterpolation (const Image& image, const glm::vec2& texCoord,  const Features& features)
{ 
    glm::vec2 texelPos { (image.width - 1) * texCoord.x, (image.height - 1) * texCoord.y};

    int clamp_texel_x = (int) texelPos.x;
    int clamp_texel_y = (int) texelPos.y;

    glm::vec2 upper_left {clamp_texel_x, clamp_texel_y}; 
    clampTexel(upper_left, image);

    glm::vec2 upper_right {upper_left.x + 1, upper_left.y};
    clampTexel(upper_right, image);

    glm::vec2 lower_left {upper_left.x, upper_left.y + 1}; 
    clampTexel(lower_left, image);

    glm::vec2 lower_right {upper_left.x + 1, upper_left.y + 1}; 
    clampTexel(lower_right, image);


    glm::vec3 upper_left_color = image.pixels[(image.height - upper_left.y - 1) * image.width + upper_left.x];
    glm::vec3 upper_right_color = image.pixels[(image.height - upper_right.y - 1) * image.width + upper_right.x];
    glm::vec3 lower_left_color = image.pixels[(image.height - lower_left.y - 1) * image.width + lower_left.x];
    glm::vec3 lower_right_color = image.pixels[(image.height - lower_right.y - 1) * image.width + lower_right.x];


    float alpha = texelPos.x - upper_left.x;
    float betta = texelPos.y - upper_left.y;  
    
    glm::vec3 interpolated_texCoord = upper_left_color * glm::vec3{(1-alpha) * (1-betta)}
     + lower_left_color  * glm::vec3{ (1 - alpha) * betta} 
     + upper_right_color * glm::vec3{ alpha * (1 - betta) } 
     + lower_right_color * glm::vec3{alpha * betta}; 

    return interpolated_texCoord; 
}