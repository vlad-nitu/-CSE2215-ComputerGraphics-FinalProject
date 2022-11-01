#include "environment_map.h"
#include <render.h>
#include <texture.h>
#include <framework/image.h>

const Image environment = Image()

glm::vec3 getEnvironmentColor(const glm::vec3& lightDirection)
{
	return glm::vec3{ 0 };
}