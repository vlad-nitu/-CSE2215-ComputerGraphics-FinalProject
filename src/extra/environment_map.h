#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <framework/ray.h>
#include <framework/image.h>

struct Features;

extern Image environment;

glm::vec3 getEnvironmentColor(const glm::vec3& lightDirection, const Features& features);