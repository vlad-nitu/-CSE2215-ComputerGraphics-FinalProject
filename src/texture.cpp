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

glm::vec3 bilinearInterpolation (const Image& image, const glm::vec2& texCoord,  const Features& features)
{ 
    glm::vec2 texelPos { (int)((image.width - 1) * texCoord.x), (int)((image.height - 1) * texCoord.y)};

    int clamp_texel_x = texelPos.x;
    int clamp_texel_y = texelPos.y;

    glm::vec2 upper_left {clamp_texel_x, clamp_texel_y}; 
    if (upper_left.x == image.width)
        upper_left.x --;
    if (upper_left.y == image.height)
        upper_left.y --;

    glm::vec2 upper_right {upper_left.x - 1, upper_left.y};
     if (upper_right.x == image.width)
        upper_right.x --;
     if (upper_right.y == image.height)
        upper_right.y --;

    glm::vec2 lower_left {upper_left.x, upper_left.y - 1}; 
    if (lower_left.x == image.width)
        lower_left.x --;
     if (lower_left.y == image.height)
        lower_left.y --;

    glm::vec2 lower_right {upper_left.x - 1, upper_left.y - 1}; 
    if (lower_right.x == image.width)
        lower_right.x --;
     if (lower_right.y == image.height)
        lower_right.y --;

    glm::vec3 upper_left_color = image.pixels[upper_left.y * image.width + upper_left.x];
    glm::vec3 upper_right_color = image.pixels[upper_right.y * image.width + upper_right.x];
    glm::vec3 lower_left_color = image.pixels[lower_left.y * image.width + lower_left.x];
    glm::vec3 lower_right_color = image.pixels[lower_right.y * image.width + lower_right.x];


    float alpha = texelPos.x - lower_left.x;
    float betta = texelPos.y - lower_left.y;  
    
    glm::vec3 interpolated_texCoord = upper_left_color * glm::vec3{(1-alpha) * (1-betta)}
     + upper_right_color  * glm::vec3{alpha * (1 - betta)} 
     + lower_left_color * glm::vec3{(1-alpha) * betta} 
     + lower_right_color * glm::vec3{alpha * betta}; 

        return interpolated_texCoord; 
}