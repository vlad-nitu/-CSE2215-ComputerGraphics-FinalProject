#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    glm::vec3 a = v0;
    glm::vec3 b = v1;
    glm::vec3 c = v2;

    a -= p;
    b -= p;
    c -= p;

    // Compute normals for triangles: u = normal of triangle PBC; v = normal of triangle PCA; w = normal of PAB
    glm::vec3 u = glm::cross(b, c);
    glm::vec3 v = glm::cross(c, a);
    glm::vec3 w = glm::cross(a, b);

    if (glm::dot(u, v) < 0.0f || glm::dot(u, w) < 0.0f || glm::dot(v, w) < 0.0f)
        return false;
    else
        return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    if (glm::dot(plane.normal, ray.direction) == 0.0f)
        return false;
    else {
        // check notes for calculations
        float possible_t = -(plane.D + glm::dot(plane.normal, ray.origin)) * 1.0f / glm::dot(plane.normal, ray.direction);

        if (possible_t <= 0)
            return false;

        if (possible_t < ray.t) {
            ray.t = possible_t;
            return true;
        } else
            return false;
    }
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 l1 = v1 - v0;
    glm::vec3 l2 = v2 - v0;
    glm::vec3 N;
    if (glm::length(glm::cross(l1, l2)) == 0) // Degenerate triangle
        N = glm::vec3 {};
    else
        N = glm::normalize(glm::cross(l1, l2));
    // v0 point in plane -> verifies normal form, dot(N, v0)
    plane.D = -N.x * v0.x - N.y * v0.y - N.z * v0.z;
    plane.normal = N;
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(v0, v1, v2);
    glm::vec3 N = plane.normal;
    float prev_ray = ray.t;

    if (intersectRayWithPlane(plane, ray)) {
        glm::vec3 P = ray.origin + ray.t * ray.direction;
        if (pointInTriangle(v0, v1, v2, N, P))
            return true;
        else {
            ray.t = prev_ray;
            return false;
        }
    } else
        return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    glm::vec3 u = glm::normalize(ray.direction);
    glm::vec3 o = ray.origin;
    glm::vec3 c = sphere.center;
    float r = sphere.radius;

    float delta = glm::dot(u, o - c) * glm::dot(u, o - c) - (glm::length(o - c) * glm::length(o - c) - r * r);

    if (delta < 0)
        return false;
    else {
        float d1 = -glm::dot(u, o - c) + sqrt(delta);
        float d2 = -glm::dot(u, o - c) - sqrt(delta);

        if (d1 > 0 && d2 > 0) {
            if (d1 <= d2 && d1 < ray.t) {
                ray.t = d1;
                return true;
            } else if (d2 <= d1 && d2 < ray.t) {
                ray.t = d2;
                return true;
            }
        } else if (d1 > 0 && d1 < ray.t) {
            ray.t = d1;
            return true;
        } else if (d2 > 0 && d2 < ray.t) {
            ray.t = d2;
            return true;
        }

        return false;
    }
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float txmin, txmax, tymin, tymax, tzmin, tzmax;
    // if x, y or z ray coordinate is 0, then place tmin = +INF and tmax = -INF in order to not be
    // considered on LoCs 199 and 200 (we will take max from -INF and min from +INF, so we don't take these into consideration)

    if (ray.direction.x != 0) {
        txmin = (box.lower.x - ray.origin.x) / ray.direction.x;
        txmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    } else {
        txmin = std::numeric_limits<float>::max();
        txmax = std::numeric_limits<float>::min();
    }

    if (ray.direction.y != 0) {
        tymin = (box.lower.y - ray.origin.y) / ray.direction.y;
        tymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    } else {
        tymin = std::numeric_limits<float>::max();
        tymax = std::numeric_limits<float>::min();
    }

    if (ray.direction.z != 0) {
        tzmin = (box.lower.z - ray.origin.z) / ray.direction.z;
        tzmax = (box.upper.z - ray.origin.z) / ray.direction.z;
    } else {
        tzmin = std::numeric_limits<float>::max();
        tzmax = std::numeric_limits<float>::min();
    }

    float t_in_x = (txmin < txmax) ? txmin : txmax;
    float t_out_x = (txmin > txmax) ? txmin : txmax;

    float t_in_y = (tymin < tymax) ? tymin : tymax;
    float t_out_y = (tymin > tymax) ? tymin : tymax;

    float t_in_z = (tzmin < tzmax) ? tzmin : tzmax;
    float t_out_z = (tzmin > tzmax) ? tzmin : tzmax;

    float t_in = std::max(std::max(t_in_x, t_in_y), t_in_z);
    float t_out = std::min(std::min(t_out_x, t_out_y), t_out_z);

    if (t_in > t_out || t_out < 0)
        return false;
    else if (t_in < 0) // inside box
    {
        if (t_out < ray.t) {
            ray.t = t_out;
            return true;
        } else
            return false;
    } else {
        if (t_in < ray.t) {
            ray.t = t_in;
            return true;
        } else
            return false;
    }

}
