#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    int col = texCoord.x * image.width; // Convert to int so values are rounded down representing the line
    int row = texCoord.y * image.height; // Convert to int so values are rounded down representing the column

    if (col == image.width)
        col--;
    if (row == image.height)
        row--;

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


    glm::vec3 upper_left_color = image.pixels[upper_left.y * image.width + upper_left.x];
    glm::vec3 upper_right_color = image.pixels[upper_right.y * image.width + upper_right.x];
    glm::vec3 lower_left_color = image.pixels[lower_left.y * image.width + lower_left.x];
    glm::vec3 lower_right_color = image.pixels[lower_right.y * image.width + lower_right.x];


    float alpha = texelPos.x - upper_left.x;
    float betta = texelPos.y - upper_left.y;  
    
    glm::vec3 interpolated_texCoord = upper_left_color * glm::vec3{(1-alpha) * (1-betta)}
     + lower_left_color  * glm::vec3{ (1 - alpha) * betta} 
     + upper_right_color * glm::vec3{ alpha * (1 - betta) } 
     + lower_right_color * glm::vec3{alpha * betta}; 

    return interpolated_texCoord; 
}