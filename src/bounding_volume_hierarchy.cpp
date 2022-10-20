#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>

/*
* Given a node, computer the bounding volume
*/
void BoundingVolumeHierarchy::computeAABB(Node& node)
{
    glm::vec3 low = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 high = glm::vec3 { -std::numeric_limits<float>::max() };

    for (int i = 0; i< node.children.size(); i++) {
        int primitiveIndex = node.children[i];

        // Find the mesh which contains this triangle
        int mesh = 0;
        while (m_pScene->meshes[mesh].triangles.size() <= primitiveIndex) {
            primitiveIndex -= m_pScene->meshes[mesh++].triangles.size();
        }

        // Check if the primitive is a sphere or a triangle
        if (m_pScene->meshes.size() <= mesh) { // !!!!!!!!!!!!!!!!! Make sure the equality sign is right
            Sphere sphere = m_pScene->spheres[primitiveIndex];

            // Computer the maximum and minimum coordinates of the sphere
            glm::vec3 sphereMin = sphere.center - glm::vec3 { sphere.radius };
            glm::vec3 sphereMax = sphere.center + glm::vec3 { sphere.radius };

            // Update the max and min values coordinate-wise
            low = glm::min(sphereMin, low);
            high = glm::max(sphereMax, high);
        } else {
            // Get the triangle from the coresponding mesh
            glm::uvec3 tri = m_pScene->meshes[mesh].triangles[primitiveIndex];

            // Get the vertices coordinates of the triangle
            glm::vec3 v0 = m_pScene->meshes[mesh].vertices[tri[0]].position;
            glm::vec3 v1 = m_pScene->meshes[mesh].vertices[tri[1]].position;
            glm::vec3 v2 = m_pScene->meshes[mesh].vertices[tri[2]].position;

            // Update the max and min values coordinate-wise
            low = glm::min(v0, low);
            high = glm::max(v0, high);

            low = glm::min(v1, low);
            high = glm::max(v1, high);

            low = glm::min(v2, low);
            high = glm::max(v2, high);
        }
    }

    node.lower = low;
    node.upper = high;
}

/// <summary>
/// Computes the AABB for a given node and decides to split it if the node has more than one primitive inside
/// </summary>
/// <param name="node"> The node that is curently being subdivided. </param>
/// <param name="centroids"> The list containing the centroids for all the primitives. </param>
/// <param name="axis"> The axis on which the node will be split. 0 for x, 1 for y, 2 for z. </param>
/// <param name="depth"> The depth of the given node. </param>
void BoundingVolumeHierarchy::subdivideNode(Node& node, std::vector<glm::vec3> centroids, int axis, int depth)
{
    computeAABB(node);

    // Check if we have reached a new bigger depth in the tree (levels = max_depth + 1 since root has depth = 0)
    m_numLevels = std::max(m_numLevels, depth + 1);

    if (node.children.size() <= 5) {
        nodes.push_back(node); // Add the node to the hierarchy

        return;
    } else {
        // We further need to subdivide
        node.isLeaf = false;
        m_numLeaves--; // Current node is no longer a leaf

        Node leftChild = Node(true);
        Node rightChild = Node(true);
        m_numLeaves += 2; // Children are by default leafs

        // Sort method taken from https://en.cppreference.com/w/cpp/algorithm/sort
        /*
        * Sort the indices based on the centroids by the current axis in order to find the median one.
        */
        std::sort(node.children.begin(), node.children.end(), [centroids, axis](int a, int b) {
            if (axis == 0)
                return centroids[a].x < centroids[b].x;
            else if (axis == 1)
                return centroids[a].y < centroids[b].y;
            else
                return centroids[a].z < centroids[b].z;
        });

        // Move the primitives to the children based on the centroids
            for (int i = 0; i < node.children.size(); i++) {

                if (i < node.children.size() / 2) {
                    leftChild.children.push_back(node.children[i]);
                } else {
                    /*if (axis == 0 && centroids[node.children[i]].x == centroids[node.children[node.children.size() / 2 - 1]].x) {
                        leftChild.children.push_back(node.children[i]);
                    } else if (axis == 1 && centroids[node.children[i]].y == centroids[node.children[node.children.size() / 2 - 1]].y) {
                        leftChild.children.push_back(node.children[i]);
                    } else if (axis == 2 && centroids[node.children[i]].z == centroids[node.children[node.children.size() / 2 - 1]].z) {
                        leftChild.children.push_back(node.children[i]);
                    } else {*/
                        rightChild.children.push_back(node.children[i]);
                    //}
                }
            }
        node.children.clear(); // remove the children index from the parent

        subdivideNode(leftChild, centroids, (axis + 1) % 3, depth + 1);
        node.children.push_back(nodes.size() - 1); // Remember the index of the left child for the parent
        
        subdivideNode(rightChild, centroids, (axis + 1) % 3, depth + 1);
        node.children.push_back(nodes.size() - 1); // remember the index of the right child for the parent
        
        // Add the node to the hierarchy
        nodes.push_back(node);
    }
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    m_numLeaves = 1;
    m_numLevels = 1;

    /*
    * Create a list containing the centers of all the triangles and spheres
    * 
    * Follows the same ordering as for the flattened list of triangles used in the node struct
    */
    std::vector<glm::vec3> centroids;

    /*
     * PLS pay attention to this
     */
    Node root = Node(true); // No idea if I need to use the 'new' keyword or not

    int primitiveIndex = 0;
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            const auto v0 = mesh.vertices[tri[0]].position;
            const auto v1 = mesh.vertices[tri[1]].position;
            const auto v2 = mesh.vertices[tri[2]].position;

            glm::vec3 centre = (v0 + v1 + v2) / glm::vec3 {3};
            centroids.push_back(centre);

            // Add all the primitives inside of the node
            root.children.push_back(primitiveIndex++);
        }
    }

    for (const auto& sphere : m_pScene->spheres) {
        centroids.push_back(sphere.center);

        // Add the sphere primitives inside the node
        root.children.push_back(primitiveIndex++);
    }

    subdivideNode(root, centroids, 0, 0);
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

void BoundingVolumeHierarchy::showLevel(Node& node, int currentLevel, int targetLevel, std::vector <AxisAlignedBox>& toDraw)
{
    if (currentLevel == targetLevel) {
        AxisAlignedBox aabb {node.lower, node.upper};

        toDraw.push_back(aabb);
    } else {
        if (!node.isLeaf) {
            showLevel(nodes[node.children[0]], currentLevel + 1, targetLevel, toDraw);
            showLevel(nodes[node.children[1]], currentLevel + 1, targetLevel, toDraw);
        }
    }
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    std::vector<AxisAlignedBox> toDraw;
    showLevel(nodes[nodes.size() - 1], 0, level, toDraw);

    for (const auto& aabb : toDraw) {
        drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1), 0.4f);
    }
}

void BoundingVolumeHierarchy::getLeaf(int index, int& leafIdx, int& result)
{
    Node node = nodes[index];

    if (node.isLeaf) {
        if (leafIdx == 0)
            result = index;

        leafIdx--;
    } else {
        if (!node.isLeaf) {
            getLeaf(node.children[0], leafIdx, result);
            getLeaf(node.children[1], leafIdx, result);
        }
    }
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    int leafIndex = nodes.size() - 1;
    int result = 0;

    getLeaf(leafIndex, leafIdx, result);

    drawAABB(AxisAlignedBox { nodes[result].lower, nodes[result].upper }, DrawMode::Wireframe, glm::vec3(1), 0.4f);

    for (int i = 0; i < nodes[result].children.size(); i++) {
        int primitiveIndex = nodes[result].children[i];

        // Find the mesh which contains this triangle
        int mesh = 0;
        while (m_pScene->meshes[mesh].triangles.size() <= primitiveIndex) {
            primitiveIndex -= m_pScene->meshes[mesh++].triangles.size();
        }

        // Check if the primitive is a sphere or a triangle
        if (m_pScene->meshes.size() <= mesh) {
            Sphere sphere = m_pScene->spheres[primitiveIndex];

            drawSphere(sphere);
        } else {
            // Get the triangle from the coresponding mesh
            glm::uvec3 tri = m_pScene->meshes[mesh].triangles[primitiveIndex];

            // Get the vertices coordinates of the triangle
            auto v0 = m_pScene->meshes[mesh].vertices[tri[0]];
            auto v1 = m_pScene->meshes[mesh].vertices[tri[1]];
            auto v2 = m_pScene->meshes[mesh].vertices[tri[2]];

            drawTriangle(v0, v1, v2, glm::vec3 { 0.3f, 0.4f, 0.8f });
        }
    }

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
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}