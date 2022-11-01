#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include <glm/glm.hpp>
#include <queue>

// Draw the debug for the interpolated normals
bool drawNormalInterpolationDebug = false;

// Draw the debug for BVH traversal
bool rayNodeIntersectionDebug = false;
bool drawUnvisited = false;
bool drawSAH_Debug = false;

/// <summary>
/// Given a primitive's index it will update the minimum and maximum coordinates of an AABBs corners
/// </summary>
/// <param name="primitiveIndex"> The index of the primitive that is used to update</param>
/// <param name="low"> The lower boundary of the AABB </param>
/// <param name="high"> The upper boundary of the AABB </param>
void BoundingVolumeHierarchy::updateAABB(int primitiveIndex, glm::vec3& low, glm::vec3& high)
{
    // Find the mesh which contains this triangle
    int mesh = 0;
    while (m_pScene->meshes.size() > 0 && m_pScene->meshes[mesh].triangles.size() <= primitiveIndex) {
        primitiveIndex -= m_pScene->meshes[mesh++].triangles.size();
    }

    // Check if the primitive is a sphere or a triangle
    if (m_pScene->meshes.size() <= mesh) {
        Sphere& sphere = m_pScene->spheres[primitiveIndex];

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
        glm::vec3& v0 = m_pScene->meshes[mesh].vertices[tri[0]].position;
        glm::vec3& v1 = m_pScene->meshes[mesh].vertices[tri[1]].position;
        glm::vec3& v2 = m_pScene->meshes[mesh].vertices[tri[2]].position;

        // Update the max and min values coordinate-wise
        low = glm::min(low, glm::min(v0, glm::min(v1, v2)));
        high = glm::max(high, glm::max(v0, glm::max(v1, v2)));
    }
}

/// <summary>
/// Computes the AABB for a given node and decides to split it if the node has more than one primitive inside
/// </summary>
/// <param name="node"> The node that is curently being subdivided. </param>
/// <param name="centroids"> The list containing the centroids for all the primitives. </param>
/// <param name="axis"> The axis on which the node will be split. 0 for x, 1 for y, 2 for z. </param>
/// <param name="depth"> The depth of the given node. </param>
void BoundingVolumeHierarchy::subdivideNode(Node& node, const std::vector<glm::vec3>& centroids, int axis, int depth)
{
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
        std::sort(node.children.begin(), node.children.end(), [&centroids, axis](int a, int b) {
            if (axis == 0)
                return centroids[a].x < centroids[b].x;
            else if (axis == 1)
                return centroids[a].y < centroids[b].y;
            else
                return centroids[a].z < centroids[b].z;
        });

        // Move the primitives to the children based on the centroids
        glm::vec3 leftMin = glm::vec3 { std::numeric_limits<float>::max() };
        glm::vec3 leftMax = glm::vec3 { -std::numeric_limits<float>::max() };

        glm::vec3 rightMin = glm::vec3 { std::numeric_limits<float>::max() };
        glm::vec3 rightMax = glm::vec3 { -std::numeric_limits<float>::max() };

        for (int i = 0; i < node.children.size(); i++) {

            if (i < node.children.size() / 2) {
                leftChild.children.push_back(node.children[i]);
                updateAABB(node.children[i], leftMin, leftMax); // Update the AABB of the left child
            } else {
                rightChild.children.push_back(node.children[i]);
                updateAABB(node.children[i], rightMin, rightMax); // Update the AABB of the right child
            }
        }

        // Assign AABBs to children nodes
        leftChild.lower = leftMin;
        leftChild.upper = leftMax;

        rightChild.lower = rightMin;
        rightChild.upper = rightMax;

        node.children.clear(); // Remove the children index from the parent

        subdivideNode(leftChild, centroids, (axis + 1) % 3, depth + 1);
        node.children.push_back(nodes.size() - 1); // Remember the index of the left child for the parent

        subdivideNode(rightChild, centroids, (axis + 1) % 3, depth + 1);
        node.children.push_back(nodes.size() - 1); // remember the index of the right child for the parent

        // Add the node to the hierarchy
        nodes.push_back(node);
    }
}

void BoundingVolumeHierarchy::subdivideNodeSah(Node& node, const std::vector<AxisAlignedBox>& AABBs, const std::vector<glm::vec3>& centroids, int depth)
{
    // Check if we have reached a new bigger depth in the tree (levels = max_depth + 1 since root has depth = 0)
    m_numLevels = std::max(m_numLevels, depth + 1);

    if (node.children.size() <= MAX_PRIMITIVES_PER_LEAF || depth >= MAX_DEPTH) {
        nodes.push_back(node); // Add the node to the hierarchy
        return;
    }

    // compute AABB of centroids of triangles
    glm::vec3 low_aabb = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 high_aabb = glm::vec3 { -std::numeric_limits<float>::max() };

    int primitives = node.children.size();

    for (int i = 0; i < primitives; ++i) {
        updateAABB_SAH(centroids[node.children[i]], low_aabb, high_aabb);
    }

    float max_val = -std::numeric_limits<float>::max();
    int axis = -1;
    float cbmin, cbmax;

    if (high_aabb.x - low_aabb.x > max_val) {
        max_val = high_aabb.x - low_aabb.x;
        cbmin = low_aabb.x;
        cbmax = high_aabb.x;
        axis = 0;
    }

    if (high_aabb.y - low_aabb.y > max_val) {
        max_val = high_aabb.y - low_aabb.y;
        cbmin = low_aabb.y;
        cbmax = high_aabb.y;
        axis = 1;
    }

    if (high_aabb.z - low_aabb.z > max_val) {
        max_val = high_aabb.z - low_aabb.z;
        cbmin = low_aabb.z;
        cbmax = high_aabb.z;
        axis = 2;
    }

    std::vector<SahAABB> bins(N_BINS);
    float centr;

    // Place every triangle in the corresponding bin
    for (int i = 0; i < primitives; ++i) {
        if (axis == 0)
            centr = centroids[node.children[i]].x;
        else if (axis == 1)
            centr = centroids[node.children[i]].y;
        else
            centr = centroids[node.children[i]].z;

        int bin_index = (1 - EPSILON) * N_BINS * (centr - cbmin) / (cbmax - cbmin);
        updateAABB(node.children[i], bins[bin_index].bounds.lower, bins[bin_index].bounds.upper);
        bins[bin_index].primitiveIndexes.push_back(i);
    }

    // Store volume and # of primitives for each bin
    std::pair<float, int> defaulted { 0.0f, 0 };
    std::vector<std::pair<float, int>> costs(N_BINS, defaulted);

    // Precompute left->right then use precomputation when traversing right->left in order to obtain linear time complexity
    costs[0] = (std::make_pair(volume(bins[0].bounds), bins[0].primitiveIndexes.size()));
    AxisAlignedBox updated_box = bins[0].bounds;
    int added_triangles = bins[0].primitiveIndexes.size();

    for (int idx = 1; idx < N_BINS; ++idx) { // left->right

        if (!bins[idx].primitiveIndexes.empty()) // compute cost
        {
            unionBoxes(updated_box, bins[idx].bounds);
            added_triangles += bins[idx].primitiveIndexes.size();

            std::pair<float, int> curr_cost = { volume(updated_box), added_triangles };
            costs[idx] = curr_cost;
        }
    }

    int split_idx = -1;
    float min_val = std::numeric_limits<float>::max();

    updated_box = bins[N_BINS - 1].bounds;
    added_triangles = bins[N_BINS - 1].primitiveIndexes.size();

    for (int idx = N_BINS - 2; idx >= 0; --idx) {

        if (!bins[idx].primitiveIndexes.empty()) // compute cost if there exist primitives inside current bin
        {
            unionBoxes(updated_box, bins[idx].bounds);
            added_triangles += bins[idx].primitiveIndexes.size();

            float vol = volume(updated_box);

            float cost = costs[idx].first * costs[idx].second + vol * added_triangles;

            if (costs[idx].first * costs[idx].second <= OFFSET || vol * added_triangles <= OFFSET) // if there are no primitives to the left OR to the right of the bin -> do not perform split
                continue;

            if (cost < min_val) // both splits are not empty
            {
                min_val = cost;
                split_idx = idx;
            }
        }
    }

    if (split_idx == -1 || min_val == std::numeric_limits<float>::max()) {
        node.isLeaf = true;
        return;
    }

    // We further need to subdivide
    node.isLeaf = false;
    m_numLeaves--; // Current node is no longer a leaf

    Node leftChild = Node(true);
    Node rightChild = Node(true);
    m_numLeaves += 2; // Children are by default leafs

    // Move the primitives to the children based on the centroids
    glm::vec3 leftMin = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 leftMax = glm::vec3 { -std::numeric_limits<float>::max() };

    glm::vec3 rightMin = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 rightMax = glm::vec3 { -std::numeric_limits<float>::max() };

    for (int i = 0; i <= split_idx; i++)
        for (int j = 0; j < bins[i].primitiveIndexes.size(); ++j) {
            leftChild.children.push_back(node.children[bins[i].primitiveIndexes[j]]);
            updateAABB(node.children[bins[i].primitiveIndexes[j]], leftMin, leftMax); // Update the AABB of the left child
        }

    for (int i = split_idx + 1; i < N_BINS; ++i)
        for (int j = 0; j < bins[i].primitiveIndexes.size(); ++j) {
            rightChild.children.push_back(node.children[bins[i].primitiveIndexes[j]]);
            updateAABB(node.children[bins[i].primitiveIndexes[j]], rightMin, rightMax); // Update the AABB of the left child
        }

    // Assign AABBs to children nodes
    leftChild.lower = leftMin;
    leftChild.upper = leftMax;

    rightChild.lower = rightMin;
    rightChild.upper = rightMax;

    node.children.clear(); // Remove the children index from the parent

    subdivideNodeSah(leftChild, AABBs, centroids, depth + 1);
    node.children.push_back(nodes.size() - 1); // Remember the index of the left child for the parent

    subdivideNodeSah(rightChild, AABBs, centroids, depth + 1);
    node.children.push_back(nodes.size() - 1); // Remember the index of the right child for the parent

    nodes.push_back(node);
}

void BoundingVolumeHierarchy::updateAABB_SAH(const glm::vec3& v, glm::vec3& lower, glm::vec3& upper)
{
    lower = glm::min(v, lower);
    upper = glm::max(v, upper);
    return;
}

void BoundingVolumeHierarchy::unionBoxes(AxisAlignedBox& updated_box, const AxisAlignedBox& next_box)
{
    updated_box.lower = glm::min(updated_box.lower, glm::min(next_box.lower, next_box.upper));
    updated_box.upper = glm::max(updated_box.upper, glm::max(next_box.lower, next_box.upper));
    return;
}

float BoundingVolumeHierarchy::volume(const AxisAlignedBox& AABB)
{
    glm::vec3 diff = AABB.upper - AABB.lower;
    return (diff.x * diff.y * diff.z);
}

glm::vec3 BoundingVolumeHierarchy::computeAABB_centroid(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    glm::vec3 lower = glm::min(v0, glm::min(v1, v2));
    glm::vec3 upper = glm::max(v0, glm::max(v1, v2));
    AxisAlignedBox aabb = { lower, upper };
    return (lower + upper) / glm::vec3 { 2 };
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{
    m_numLeaves = 1;
    m_numLevels = 1;

    /*
     * Create a list containing the centers of all the triangles and spheres
     *
     * Follows the same ordering as for the flattened list of triangles used in the node struct
     */
    std::vector<glm::vec3> centroids {};

    // Stores an AABB per primitive after encapsulating each into an AABB
    std::vector<AxisAlignedBox> AABBs {};

    Node root = Node(true);

    glm::vec3 low = glm::vec3 { std::numeric_limits<float>::max() };
    glm::vec3 high = glm::vec3 { -std::numeric_limits<float>::max() };

    int primitiveIndex = 0;
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            const auto& v0 = mesh.vertices[tri[0]].position;
            const auto& v1 = mesh.vertices[tri[1]].position;
            const auto& v2 = mesh.vertices[tri[2]].position;

            if (features.extra.enableBvhSahBinning) {
                // Compute AABB for each triangle
                // centroids.push_back(computeAABB_centroid(v0, v1, v2));
                // glm::vec3 centre = (v0 + v1 + v2) / glm::vec3 { 3 };

                centroids.push_back(computeAABB_centroid(v0, v1, v2)); // centroid computed on AABB centroid instead of triangle's centroid

                glm::vec3 low_aabb = glm::vec3 { std::numeric_limits<float>::max() };
                glm::vec3 high_aabb = glm::vec3 { -std::numeric_limits<float>::max() };

                updateAABB(primitiveIndex, low_aabb, high_aabb);

                AABBs.push_back({ low_aabb, high_aabb });
            } else {
                // Compute centroids
                glm::vec3 centre = (v0 + v1 + v2) / glm::vec3 { 3 };
                centroids.push_back(centre);
            }

            // Compute root's AABB boundry
            low = glm::min(low, glm::min(v0, glm::min(v1, v2)));
            high = glm::max(high, glm::max(v0, glm::max(v1, v2)));

            // Add all the primitives inside of the node
            root.children.push_back(primitiveIndex++);
        }
    }

    for (const auto& sphere : m_pScene->spheres) {
        centroids.push_back(sphere.center);

        // Computer the maximum and minimum coordinates of the sphere
        glm::vec3 sphereMin = sphere.center - glm::vec3 { sphere.radius };
        glm::vec3 sphereMax = sphere.center + glm::vec3 { sphere.radius };

        // Update the max and min values coordinate-wise
        low = glm::min(sphereMin, low);
        high = glm::max(sphereMax, high);

        // Add the sphere primitives inside the node
        root.children.push_back(primitiveIndex++);
    }

    root.lower = low;
    root.upper = high;

    if (features.extra.enableBvhSahBinning)
        subdivideNodeSah(root, AABBs, centroids, 0);
    else
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

void BoundingVolumeHierarchy::showLevel(const Node& node, int currentLevel, int targetLevel)
{
    if (currentLevel == targetLevel) {
        AxisAlignedBox aabb { node.lower, node.upper };

        drawAABB(aabb, DrawMode::Wireframe, glm::vec3(1), 0.4f);
    } else {
        if (!node.isLeaf) {
            showLevel(nodes[node.children[0]], currentLevel + 1, targetLevel);
            showLevel(nodes[node.children[1]], currentLevel + 1, targetLevel);
        }
    }
}

void BoundingVolumeHierarchy::showLevelSAH(const Node& node, int currentLevel, int targetLevel)
{
    if (currentLevel == targetLevel) {

        if (!node.isLeaf) // make sure that children exist
        {
            AxisAlignedBox aabb_left { nodes[node.children[0]].lower, nodes[node.children[0]].upper };
            AxisAlignedBox aabb_right { nodes[node.children[1]].lower, nodes[node.children[1]].upper };

            drawAABB(aabb_left, DrawMode::Wireframe, glm::vec3 { 0, 1, 0 }, 0.4f);
            drawAABB(aabb_right, DrawMode::Wireframe, glm::vec3 { 1, 0, 0 }, 0.4f);
        }

    } else if (!node.isLeaf) {
        // only recurse on left child for debug purpose
        showLevelSAH(nodes[node.children[0]], currentLevel + 1, targetLevel);
    }
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{

    if (!drawSAH_Debug) // Classic BVH
        showLevel(nodes[nodes.size() - 1], 0, level);
    else // SAH
        showLevelSAH(nodes[nodes.size() - 1], 0, level);
}

void BoundingVolumeHierarchy::getLeaf(int index, int& leafIdx, int& result)
{
    Node& node = nodes[index];

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
        while (m_pScene->meshes.size() > 0 && m_pScene->meshes[mesh].triangles.size() <= primitiveIndex) {
            primitiveIndex -= m_pScene->meshes[mesh++].triangles.size();
        }

        // Check if the primitive is a sphere or a triangle
        if (m_pScene->meshes.size() <= mesh) {
            Sphere& sphere = m_pScene->spheres[primitiveIndex];

            drawSphere(sphere);
        } else {
            // Get the triangle from the coresponding mesh
            glm::uvec3 tri = m_pScene->meshes[mesh].triangles[primitiveIndex];

            // Get the vertices coordinates of the triangle
            auto& v0 = m_pScene->meshes[mesh].vertices[tri[0]];
            auto& v1 = m_pScene->meshes[mesh].vertices[tri[1]];
            auto& v2 = m_pScene->meshes[mesh].vertices[tri[2]];

            drawTriangle(v0, v1, v2, glm::vec3 { 1.0f, 0.2f, 0.6f });
        }
    }

    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}

/// <summary>
/// Given a primitive's index it will draw it in a different color or it's interpolated normals
/// </summary>
/// <param name="primitiveIndex"> The index of the primitive that is to be drawn </param>
void BoundingVolumeHierarchy::drawPrimitive(int primitiveIndex, const Ray& ray, const Features& features, const HitInfo& hitInfo) const
{
    // Find the mesh which contains this triangle
    int mesh = 0;
    while (m_pScene->meshes.size() > 0 && m_pScene->meshes[mesh].triangles.size() <= primitiveIndex) {
        primitiveIndex -= m_pScene->meshes[mesh++].triangles.size();
    }

    // Check if the primitive is a sphere or a triangle
    if (m_pScene->meshes.size() <= mesh) {
        Sphere& sphere = m_pScene->spheres[primitiveIndex];

        if (rayNodeIntersectionDebug)
            drawSphere(sphere);

        if (features.enableNormalInterp && drawNormalInterpolationDebug) {
            glm::vec3 p = ray.origin + ray.t * ray.direction;

            drawRay({ p, hitInfo.normal, 0.2f }, glm::vec3 { 0, 1, 0 }); // Used normal
            drawRay({ p, p - sphere.center, 0.4f }, glm::vec3 { 0, 0, 1 }); // Actual normal normal
        }

    } else {
        // Get the triangle from the coresponding mesh
        glm::uvec3 tri = m_pScene->meshes[mesh].triangles[primitiveIndex];

        // Get the vertices coordinates of the triangle
        auto& v0 = m_pScene->meshes[mesh].vertices[tri[0]];
        auto& v1 = m_pScene->meshes[mesh].vertices[tri[1]];
        auto& v2 = m_pScene->meshes[mesh].vertices[tri[2]];

        if (rayNodeIntersectionDebug)
            drawTriangle(v0, v1, v2, glm::vec3 { 1.0f, 0.2f, 0.6f });

        if (features.enableNormalInterp && drawNormalInterpolationDebug) {
            drawRay({ v0.position, v0.normal, 0.2f }, glm::vec3 { 0, 0, 1 });
            drawRay({ v1.position, v1.normal, 0.2f }, glm::vec3 { 0, 0, 1 });
            drawRay({ v2.position, v2.normal, 0.2f }, glm::vec3 { 0, 0, 1 });

            drawRay({ ray.origin + ray.t * ray.direction, hitInfo.normal, 0.2f }, glm::vec3 { 0, 1, 0 }); // Normal
            drawRay({ ray.origin + ray.t * ray.direction, interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord), 0.4f }, glm::vec3 { 0, 0, 1 }); // Interpolated normal
        }
    }
}

/// <summary>
/// Checks if a given ray starts inside of the given AABB
/// </summary>
/// <param name="ray"> The ray that is to be checked </param>
/// <param name="aabb"> The AABB which is to be determined if it contains the origin of the ray</param>
/// <returns> True if the origin of the ray is inside the given AABB, false otherwise</returns>
bool BoundingVolumeHierarchy::isInAABB(const Ray& ray, const AxisAlignedBox& aabb) const
{
    return ((aabb.lower.x <= ray.origin.x && ray.origin.x <= aabb.upper.x) && (aabb.lower.y <= ray.origin.y && ray.origin.y <= aabb.upper.y) && (aabb.lower.z <= ray.origin.z && ray.origin.z <= aabb.upper.z));
}

/// <summary>
/// For a leaf node checks the intersection between this node's primitives and the given ray and, if an intersection happens, updates hitInfo
/// </summary>
/// <param name="node"> The node whose primitives are to be checked </param>
/// <param name="ray"> The ray with which the intersection will happen </param>
/// <param name="hitInfo"> Struct containg information about the intersection(normals, barycentric coordinates, etc) </param>
/// <param name="features"> Struct containg info about features </param>
/// <param name="bestPrimitiveIndex"> The index of the curently intersected primitive. Is to be updated if a closer to the origin of the ray one is found </param>
/// <returns> True if an intersection happens, false otherwise </returns>
bool BoundingVolumeHierarchy::testPrimitives(const Node& node, Ray& ray, HitInfo& hitInfo, const Features& features, int& bestPrimitiveIndex) const
{
    bool hit = false;

    // Storing information about the best triangle intersection
    Vertex bestV0;
    Vertex bestV1;
    Vertex bestV2;
    if (m_pScene->meshes.size() > 0) {
        bestV0 = m_pScene->meshes[0].vertices[0];
        bestV1 = m_pScene->meshes[0].vertices[0];
        bestV2 = m_pScene->meshes[0].vertices[0];
    } else {
        bestV0 = Vertex();
        bestV1 = Vertex();
        bestV2 = Vertex();
    }
    float bestTriangleT = ray.t;
    int bestTriangleIndex = -1;

    // Storing information about the best sphere intersection
    Sphere bestSphere;
    if (m_pScene->spheres.size() > 0)
        bestSphere = m_pScene->spheres[0];
    else
        bestSphere = Sphere();
    float bestSphereT = ray.t;
    int bestSphereIndex = -1;

    // Iterating over all this node's primitives
    for (int i = 0; i < node.children.size(); i++) {
        int primitiveIndex = node.children[i];

        // Find the mesh which contains this triangle
        int mesh = 0;
        while (m_pScene->meshes.size() > 0 && m_pScene->meshes[mesh].triangles.size() <= primitiveIndex) {
            primitiveIndex -= m_pScene->meshes[mesh++].triangles.size();
        }

        // Check if the primitive is a sphere or a triangle
        if (m_pScene->meshes.size() <= mesh) {
            Sphere sphere = m_pScene->spheres[primitiveIndex]; // Get the sphere

            hit |= intersectRayWithShape(sphere, ray, hitInfo); // Check for interesction

            // Update info if it represents a better intersection
            if (ray.t < bestSphereT) {
                bestSphereT = ray.t;
                bestSphere = sphere;
                bestSphereIndex = node.children[i];
            }

        } else {
            // Get the triangle from the coresponding mesh
            glm::uvec3& tri = m_pScene->meshes[mesh].triangles[primitiveIndex];

            // Get the vertices coordinates of the triangle
            auto& v0 = m_pScene->meshes[mesh].vertices[tri[0]];
            auto& v1 = m_pScene->meshes[mesh].vertices[tri[1]];
            auto& v2 = m_pScene->meshes[mesh].vertices[tri[2]];

            float oldT = ray.t;
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {

                if (ray.t < bestTriangleT) {
                    bestTriangleT = ray.t;
                    bestV0 = v0;
                    bestV1 = v1;
                    bestV2 = v2;
                    bestTriangleIndex = node.children[i];
                }

                // Compute barycentric coordinates
                hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, ray.origin + ray.t * ray.direction);

                // If feature is enabled compute texture coordinates and interpolated normal
                if (features.enableNormalInterp) {
                    hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);

                    hitInfo.normal = (glm::dot(ray.direction, hitInfo.normal) > 0) ? -hitInfo.normal : hitInfo.normal;
                } else {
                    glm::vec3 normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));

                    hitInfo.normal = (glm::dot(ray.direction, normal) > 0) ? -normal : normal;
                }

                if (features.enableTextureMapping) {
                    hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                }

                hitInfo.material = m_pScene->meshes[mesh].material;
                hit = true;
            } else {
                ray.t = oldT;
            }
        }
    }

    // Return the index of the hit primitive
    if (hit && bestSphereT < bestTriangleT && bestSphereIndex != -1) {
        if (bestSphereT <= ray.t) {
            glm::vec3 p = ray.origin + ray.t * ray.direction;
            hitInfo.normal = glm::normalize(p - bestSphere.center);

            hitInfo.normal = (glm::dot(ray.direction, hitInfo.normal) > 0) ? -hitInfo.normal : hitInfo.normal;

            hitInfo.material = bestSphere.material;
            bestPrimitiveIndex = bestSphereIndex;
        }
    } else if (bestTriangleIndex != -1) {
        if (bestTriangleT <= ray.t)
            bestPrimitiveIndex = bestTriangleIndex;
    }

    return hit;
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

        Vertex bestV0;
        Vertex bestV1;
        Vertex bestV2;
        if (m_pScene->meshes.size() > 0) {
            bestV0 = m_pScene->meshes[0].vertices[0];
            bestV1 = m_pScene->meshes[0].vertices[0];
            bestV2 = m_pScene->meshes[0].vertices[0];
        } else {
            bestV0 = Vertex();
            bestV1 = Vertex();
            bestV2 = Vertex();
        }
        float bestTriangleT = ray.t;

        Sphere bestSphere;
        if (m_pScene->spheres.size() > 0)
            bestSphere = m_pScene->spheres[0];
        else
            bestSphere = Sphere();
        float bestSphereT = ray.t;

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

                    } else {
                        glm::vec3 normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));

                        hitInfo.normal = (glm::dot(ray.direction, normal) > 0) ? -normal : normal;
                    }

                    if (features.enableTextureMapping) {
                        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
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
        if (features.enableNormalInterp && drawNormalInterpolationDebug) {
            if (bestTriangleT < bestSphereT) {
                drawRay({ bestV0.position, bestV0.normal, 0.2f }, glm::vec3 { 0, 0, 1 });
                drawRay({ bestV1.position, bestV1.normal, 0.2f }, glm::vec3 { 0, 0, 1 });
                drawRay({ bestV2.position, bestV2.normal, 0.2f }, glm::vec3 { 0, 0, 1 });

                drawRay({ ray.origin + ray.t * ray.direction, hitInfo.normal, 0.2f }, glm::vec3 { 0, 1, 0 }); // Normal
                drawRay({ ray.origin + ray.t * ray.direction, interpolateNormal(bestV0.normal, bestV1.normal, bestV2.normal, hitInfo.barycentricCoord), 0.4f }, glm::vec3 { 0, 0, 1 }); // Interpolated normal
            } else {
                glm::vec3 p = ray.origin + ray.t * ray.direction;
                hitInfo.normal = glm::normalize(p - bestSphere.center);

                hitInfo.normal = (glm::dot(ray.direction, hitInfo.normal) > 0) ? -hitInfo.normal : hitInfo.normal;
                hitInfo.material = bestSphere.material;

                drawRay({ p, hitInfo.normal, 0.2f }, glm::vec3 { 0, 1, 0 }); // Used normal
                drawRay({ p, p - bestSphere.center, 0.4f }, glm::vec3 { 0, 0, 1 }); // Actual normal normal
            }
        }

        return hit;
    } else {
        struct Trav {
            float t;
            int NodeIndex;
        };

        auto my_comp = [](const Trav& z1, const Trav& z2) { // Sort increasing by the 1st argument and decreasing by the second one
            return z1.t > z2.t;
        };

        std::priority_queue<Trav, std::vector<Trav>, decltype(my_comp)> queue(my_comp); // Create priority queue
        bool hitPrimitive = false;

        AxisAlignedBox rootAABB = { nodes[nodes.size() - 1].lower, nodes[nodes.size() - 1].upper };
        float oldT = ray.t;

        int hitIndex = -1;

        if (intersectRayWithShape(rootAABB, ray) || isInAABB(ray, rootAABB)) {
            if (isInAABB(ray, rootAABB))
                queue.push({ 0, (int)(nodes.size() - 1) }); // Push root to queue with key 0 (matters in case both origin and end are inside root AABB)
            else
                queue.push({ ray.t, (int)(nodes.size() - 1) }); // Push root to queue
            ray.t = oldT;

            // Go through the tree until no better options are available
            while (queue.size() > 0) {
                const Trav current = queue.top();
                queue.pop();
                const Node& node = nodes[current.NodeIndex];

                // Draw the AABBs that the ray intersects
                if (rayNodeIntersectionDebug)
                    drawAABB(AxisAlignedBox { node.lower, node.upper }, DrawMode::Wireframe, glm::vec3 { 1 }, 0.4f);

                // The current node has potential for a better intersection
                if (current.t <= ray.t) {

                    if (node.isLeaf) {
                        hitPrimitive |= testPrimitives(node, ray, hitInfo, features, hitIndex);
                    } else {

                        // Node has left child
                        if (node.children.size() > 0) {
                            const AxisAlignedBox leftChildBox = { nodes[node.children[0]].lower, nodes[node.children[0]].upper };
                            oldT = ray.t;

                            if (isInAABB(ray, leftChildBox)) {

                                queue.push({ 0, node.children[0] }); // Insert the node as a primitive can be intersected at any time

                            } else if (intersectRayWithShape(leftChildBox, ray)) {

                                // Insert into queue only if the node has potential for a better intersection
                                if (ray.t <= oldT)
                                    queue.push({ ray.t, node.children[0] });

                                // Reset the old intersection value
                                ray.t = oldT;
                            }
                        }

                        // Node has right child
                        if (node.children.size() > 1) {
                            const AxisAlignedBox rightChildBox = { nodes[node.children[1]].lower, nodes[node.children[1]].upper };
                            oldT = ray.t;

                            if (isInAABB(ray, rightChildBox)) {

                                queue.push({ 0, node.children[1] }); // Insert the node as a primitive can be intersected at any time

                            } else if (intersectRayWithShape(rightChildBox, ray)) {

                                // Insert into queue only if the node has potential for a better intersection
                                if (ray.t <= oldT)
                                    queue.push({ ray.t, node.children[1] });

                                // Reset the old intersection value
                                ray.t = oldT;
                            }
                        }
                    }

                } else {
                    /*
                     * The node has been intersected but since it's farther than the current best
                     * primitive intersection we have no reason to visit it.
                     */
                    if (drawUnvisited)
                        drawAABB({ node.lower, node.upper }, DrawMode::Wireframe, glm::vec3 { 1.0f, 0.45f, 0.1f }, 0.4f);
                }
            }

            // Draw debug features for the hit primitive
            if (hitPrimitive && (rayNodeIntersectionDebug || drawNormalInterpolationDebug))
                drawPrimitive(hitIndex, ray, features, hitInfo);
        }

        return hitPrimitive;
    }
}