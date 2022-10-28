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

    // glm::vec2 left_top_px  = glm::vec2{0, 0} * texCoord; 
    // glm::vec2 right_top_px = glm::vec2{image.width - 1, 0} * texCoord; 
    // glm::vec2 left_bottom_px = glm::vec2{0, image.height - 1} * texCoord; 
    // glm::vec2 right_bottom_px = glm::vec2{image.width - 1, image.height - 1} * texCoord;

    glm::vec2 texelPos { (image.width - 1) * texCoord[0], (image.height - 1) * texCoord[1]};

    int clamp_texel_x = (int)(texelPos.x);
    int clamp_texel_y = (int)(texelPos.y);

    glm::vec2 upper_left {clamp_texel_x, clamp_texel_y}; 
    glm::vec2 upper_right {upper_left.x + 1, upper_left.y};
    glm::vec2 lower_left {upper_left.x, upper_left.y + 1}; 
    glm::vec2 lower_right {upper_left.x + 1, upper_left.y + 1}; 

    // printf("%f %f\n", upper_left.x ,upper_left.y); 

    float alpha = texelPos.x - lower_left.x;
    float betta = texelPos.y - lower_left.y;  
    
    glm::vec2 interpolated_texCoord = upper_left * glm::vec2{(1-alpha) * (1-betta)}
     + upper_right  * glm::vec2{alpha * (1 - betta)} 
     + lower_left * glm::vec2{(1-alpha) * betta} 
     + lower_right * glm::vec2{alpha * betta}; 

    // glm::vec2 rescaled {interpolated_texCoord.x / (image.width - 1), interpolated_texCoord.y / (image.height - 1)};
    //   printf("%f %f\n", interpolated_texCoord.x ,interpolated_texCoord.y); 
    //    printf("%f %f\n", rescaled.x, rescaled.y); 

        int col = interpolated_texCoord.x; 
        int row = interpolated_texCoord.y; 
        if (col == image.width)
        col--;
        if (row == image.height)
        row--;

     // printf("%d\n", row * image.width + col);

     return image.pixels[row * image.width + col];
}