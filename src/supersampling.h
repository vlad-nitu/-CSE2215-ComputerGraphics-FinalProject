Ray generateSample(Screen& screen);

glm::vec3 pixelResult(const int sampleCount, const int x, const int y, const glm::ivec2& windowResolution);

void supersampling(const int sampleCount, const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features);