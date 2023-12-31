#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;

extern bool drawNormalInterpolationDebug;

extern bool rayNodeIntersectionDebug;
extern bool drawUnvisited;

extern bool drawSAH_Debug; 

struct Node {
    bool isLeaf; // false for interior node and true for leafs.

    glm::vec3 lower; // vec3 representing the lower values for each axis.
    glm::vec3 upper; // vec3 representing the upper values for each axis.

    // Default constructor
    Node() {}

    // Basic constructor
    Node(bool leaf)
    {
        isLeaf = leaf;
    }

    // Full constructor
    Node(bool leaf, glm::vec3 low, glm::vec3 high)
    {
        isLeaf = leaf;
        lower = low;
        upper = high;
    }
    
    /*
    * Contains the indices of this node's children in their respective vector.
    * 
    * There are two cases:
    * 
    * Case 1: Node is an interior node
    * The values represent the indices of the children in the vector containing 
    * all of the nodes.
    * 
    * Case 2: Node is a leaf
    * Because primitives are divided among meshes and spheres we represent them 
    * as if contained in a single continuous vector. 
    * To get the mesh and position
    * of a primitive we substract mesh sizes one by one until the curent mesh size
    * is bigger than the remaining index of the primitive.
    * After substraction all the meshes we consider that the primitive is a sphere
    */
    std::vector<int> children; 
};

struct SahAABB { 
    AxisAlignedBox bounds {glm::vec3{std::numeric_limits<float>::max()},
                           glm::vec3{-std::numeric_limits<float>::max()}}; // boundaries of AABB
    std::vector<int> primitiveIndexes; // index in node.children vector
};


class BoundingVolumeHierarchy {
private:
    
    // Number of bins chosen for SAH+Binning
    const int N_BINS = 16;
    // Offset for calculating `bin_index` in order to remain within boundaries
    const double EPSILON = 1e-4;
    // Max # of primitives / leaf allowed for BVH
    const int MAX_PRIMITIVES_PER_LEAF = 4;
    // Max depth BVH can reach
    const int MAX_DEPTH = 15;
    // Offset for comparing floating point numbers with 0
    const double OFFSET = 1e-6;

    void updateAABB(int primitiveIndex, glm::vec3& low, glm::vec3& high);

    void subdivideNode(Node& node, const std::vector<glm::vec3>& centroids, int axis, int depth);

    void showLevel(const Node& node, int currentLevel, int targetLevel);

    void showLevelSAH(const Node& node, int currentLevel, int targetLevel);


    void getLeaf(int index, int& leafIdx, int& result);

    void drawPrimitive(int primitiveIndex, const Ray& ray, const Features& features, const HitInfo& hitInfo) const;

    bool isInAABB(const Ray& ray, const AxisAlignedBox& aabb) const;

    bool testPrimitives(const Node& node, Ray& ray, HitInfo& hitInfo, const Features& features, int& bestPrimitiveIndex) const;

    // Splits each node and performs all necessary computations for SAH+Binning step for `node`, then recurses on node.left and node.right
    void subdivideNodeSah(Node& node, const std::vector<glm::vec3>& centroids, int depth);
    
    // Computes the centroid of an AABB as `lower + upper / 2`
    glm::vec3 computeAABB_centroid(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2);

    // Calculates the volume of an AABB 
    float volume(const AxisAlignedBox& AABB);

    // Unifies `updated_box` with `next_box` in order to make a bigger box that encapsulates both
    void unionBoxes(AxisAlignedBox& updated_box, const AxisAlignedBox& next_box);

    // Updates `lower` and `upper` as being the min / max compared to `v` element-wise, respectively
    void updateAABB_SAH (const glm::vec3& v, glm::vec3& lower, glm::vec3& upper);

public:
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;

private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> nodes;
   
};