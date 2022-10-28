#include "texture.h"
#include <framework/image.h>
#include <vector> 
#include <algorithm>
#include <cmath>
#include "draw.h"

bool drawMipMapDebug = false;
std::unordered_map<Image, std::vector<Image> > map; 

/*
Create MipMap w/ 5 levels, where initial_image is on level 0 (root) 
*/

void debugDrawMipMapLevel(int level, Ray& ray) { 

    glm::vec3 center{0,0,0}; 

    // glm::vec3 curr_lower {- (2 *level), - (2 * level), - (2 * level)};
    // glm::vec3 curr_upper {(2 * level), (2 * level), (2 * level)};

    //     glm::vec3 prev_lower {curr_lower -  glm::vec3{2 * sqrt(3), 2 * sqrt(3), 2 * sqrt(3)} };
    //     glm::vec3 prev_upper {curr_upper + glm::vec3{2 * sqrt(3), 2 * sqrt(3), 2 * sqrt(3)} };

    //     drawAABB({prev_lower, prev_upper}, DrawMode::Wireframe, glm::vec3{0,1,0}, 0.4f);
    //     drawAABB({curr_lower, curr_upper}, DrawMode::Filled, glm::vec3{1,0,0}, 0.4f);
        
        // drawSphere(center, 2 * level, glm::vec3{0,1,0});
        // drawSphere(center, 2 * (level + 1),glm::vec3{1,0,0});

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

    if (features.extra.enableMipmapTextureFiltering)
        {
             if (drawMipMapDebug)
                debugDrawMipMapLevel(level, ray); 
        }

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

    if (col == image.width)
        col--;
    if (row == image.height)
        row--;

    return image.pixels[row * image.width + col];
}