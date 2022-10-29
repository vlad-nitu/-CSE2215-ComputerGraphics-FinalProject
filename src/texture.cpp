#include "texture.h"
#include <framework/image.h>
#include <vector> 
#include <algorithm>
#include <cmath>
#include "draw.h"

bool drawMipMapDebug = false;
std::unordered_map<Image, std::vector<Image> > map; 
int mipmap_max_depth;

void debugDrawMipMapLevel(int level, Ray& ray) { 

        glm::vec3 RED {1.0f, 0.0f, 0.0f}; 
        glm::vec3 GREEN {0.0f, 1.0f, 0.0f}; 
        drawSphereCustom(2 * level, RED, ray); 
        drawSphereCustom(2 * (level + 1), GREEN, ray); 

} 

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features, int level, Ray& ray)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    
    Image img = image; 

    if (features.extra.enableMipmapTextureFiltering && drawMipMapDebug)
                debugDrawMipMapLevel(level, ray); 

    int col = texCoord.x * img.width; // Convert to int so values are rounded down representing the line
    int row = texCoord.y * img.height; // Convert to int so values are rounded down representing the column

    if (col == img.width)
        col--;
    if (row == img.height)
        row--;

    return img.pixels[row * img.width + col];
}

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    
    int col = texCoord.x * image.width; // Convert to int so values are rounded down representing the line
    int row = texCoord.y * image.height; // Convert to int so values are rounded down representing the column

    row++;

    return image.pixels[(image.height - row) * image.width + col];
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

glm::vec3 bilinearInterpolation (const Image& image, const glm::vec2& texCoord,  const Features& features, int level, Ray& ray) {
        
        if (features.extra.enableMipmapTextureFiltering && drawMipMapDebug)
                debugDrawMipMapLevel(level, ray); 
        return bilinearInterpolation(image, texCoord, features);        
}