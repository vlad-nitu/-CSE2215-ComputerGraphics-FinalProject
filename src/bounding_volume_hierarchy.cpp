#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>


BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // TODO: implement BVH construction.
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return 1;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
}


// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;

        Vertex bestV0 = m_pScene->meshes[0].vertices[0];
        Vertex bestV1 = m_pScene->meshes[0].vertices[0];
        Vertex bestV2 = m_pScene->meshes[0].vertices[0];
        float bestTriangleT = std::numeric_limits<float>::max();
        
        Sphere bestSphere;
        if (m_pScene->spheres.size() > 0)
            bestSphere = m_pScene->spheres[0];
        else
            bestSphere = Sphere();
        float bestSphereT = std::numeric_limits<float>::max();

        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];

                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {

                    if (ray.t < bestTriangleT) {
                        bestTriangleT = ray.t;
                        bestV0 = v0;
                        bestV1 = v1;
                        bestV2 = v2;
                    }

                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.t * ray.direction);

                    // Check if normal interpolation is turned on
                    if (features.enableNormalInterp) {
                        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord); // Interpolate normal

                        hitInfo.normal = (glm::dot(ray.direction, hitInfo.normal) > 0) ? -hitInfo.normal : hitInfo.normal;

                        // Check if textures are turned on
                        if (features.enableTextureMapping) {
                            hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord); // Calculate texture coordinates
                        }

                    } else {
                        glm::vec3 normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));

                        hitInfo.normal = (glm::dot(ray.direction, normal) > 0) ? -normal : normal;
                    }

                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }

        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres) {
            hit |= intersectRayWithShape(sphere, ray, hitInfo);

            if (ray.t < bestSphereT) {
                bestSphereT = ray.t;
                bestSphere = sphere;
            }
        }

        // Draw debug for normal interpolation for the best primitive intersection
        if (features.enableNormalInterp) {
            if (bestTriangleT < bestSphereT) {
                drawRay({bestV0.position, bestV0.normal, 0.2f}, glm::vec3 { 0, 0, 1 });
                drawRay({bestV1.position, bestV1.normal, 0.2f}, glm::vec3 { 0, 0, 1 });
                drawRay({bestV2.position, bestV2.normal, 0.2f}, glm::vec3 { 0, 0, 1 });

                drawRay({ray.origin + ray.t * ray.direction, hitInfo.normal, 0.2f}, glm::vec3 { 0, 1, 0 });
                // Use this one in case the interpolated normal should have the same direction as the vertex ones
                //drawRay({ ray.origin + ray.t * ray.direction, interpolateNormal(bestV0.normal, bestV1.normal, bestV2.normal, hitInfo.barycentricCoord), 0.2f }, glm::vec3 { 0, 1, 0 });
            } else {
                glm::vec3 p = ray.origin + ray.t * ray.direction;

                drawRay({p, p - bestSphere.center, 0.2f}, glm::vec3 { 0, 1, 0 });
            }
        }

        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}