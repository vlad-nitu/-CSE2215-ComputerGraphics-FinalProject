#include "texture.h"
#include <framework/image.h>
#include <vector> 
#include <algorithm>
#include <cmath>

int mipmap_max_depth = 0;

/*
Create MipMap w/ 5 levels, where initial_image is on level 0 (root) 
*/
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

            glm::vec3 interpolated_pixel = (prev_image.pixels[x] + prev_image.pixels[x + 1] + prev_image.pixels[x + prev_w] + prev_image.pixels[x + prev_w + 1]) / 4.0f;  
            new_pixels.push_back(interpolated_pixel);

            if (x < prev_w)
                x += 2;
            else{
                x = 0; y += 2 * prev_w;
            }
        }

        Image new_image = Image(new_w, new_h, new_pixels);
        images.push_back(new_image);
    }
    return images;
}


glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features, int level)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    
    Image img = image; 

    if (features.extra.enableMipmapTextureFiltering)
        {
            mipmap_max_depth = std::log2(img.height);
            std::vector<Image> images = createImages(image);
            if (level > mipmap_max_depth)
                level = mipmap_max_depth;
            img = images[level];
        }

    int col = texCoord.x * img.width; // Convert to int so values are rounded down representing the line
    int row = texCoord.y * img.height; // Convert to int so values are rounded down representing the column

    if (col == img.width)
        col--;
    if (row == img.height)
        row--;

    return img.pixels[row * img.width + col];
}