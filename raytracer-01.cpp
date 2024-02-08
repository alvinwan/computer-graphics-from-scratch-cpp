/*
Raycast 01
==========
Implements a simple raycast in a scene with 3 spheres.

```bash
g++ raster-01.cpp -o raster01.out -std=c++20
./raster01.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/02-basic-raytracing.html
*/
#include "bmp.h"
#include <math.h>

const float INF = std::numeric_limits<float>::infinity();
const uint8_t BACKGROUND_COLOR[3] = {255, 255, 255};

struct Sphere {
    float center[3];
    float radius;
    uint8_t color[3];

    Sphere() {}

    Sphere(const std::initializer_list<float>& c, float r, const std::initializer_list<uint8_t>& col) {
        std::copy(c.begin(), c.end(), center);
        radius = r;
        std::copy(col.begin(), col.end(), color);
    }
};

bool put_pixel(
    uint8_t data[][3],
    int32_t width,
    int32_t height,
    int32_t x,
    int32_t y,
    const uint8_t color[3]
) {
    // Translate [-W/2, W/2] into [0, W]. Ditto for H
    x = width / 2 + x;
    y = height / 2 - y - 1;

    // Checks that the pixel is in bound
    if (x < 0 || x >= width || y < 0 || y >= height) {
        std::cerr << "Error: Attempted to write out-of-bounds pixel (" << x << ", " << y << ")." << std::endl;
        return false;
    }

    // Write into data, which is a flattened array where width * height in 1d
    int32_t offset = y * width + x;
    data[offset][0] = color[0];
    data[offset][1] = color[1];
    data[offset][2] = color[2];
    return true;
}

// Compute dot product between two 3d vectors
float dot(float v1[3], float v2[3]) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Elementwise subtraction between two 3d vectors
float* subtract(float v1[3], float v2[3]) {
    float* point = new float[3];
    point[0] = v1[0] - v2[0];
    point[1] = v1[1] - v2[1];
    point[2] = v1[2] - v2[2];
    return point;
}

// Convert 2d pixel coordinates to 3d viewport coordinates.
float* canvas_to_viewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    float* point = new float[3];
    point[0] = (float) x / width;
    point[1] = (float) y / height;
    point[2] = 1;
    return point;
}

// Computes intersection of ray with spheres. Returns solutions in terms of
// line parameter t.
float* intersect_ray_with_sphere(
    float origin[3],
    float direction[3],
    Sphere sphere
) {
    float* out = new float[2];
    float* center = subtract(origin, sphere.center);

    float a = dot(direction, direction);
    float b = 2 * dot(center, direction);
    float c = dot(center, center) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant > 0) {
        out[0] = -b + sqrt(discriminant) / (2 * a);
        out[1] = -b - sqrt(discriminant) / (2 * a);
    } else {
        out[0] = INF;
        out[1] = INF;
    }
    return out;
}

// Traces a ray against the spheres in the scene
const uint8_t* trace_ray(
    float origin[3],
    float direction[3],
    float min_t,
    float max_t,
    Sphere spheres[],
    int32_t num_spheres
) {
    float closest_t = INF;
    Sphere closest_sphere;

    for (int i = 0; i < num_spheres; i++) {
        float* ts = intersect_ray_with_sphere(origin, direction, spheres[i]);
        if (ts[0] < closest_t && min_t < ts[0] && ts[0] < max_t) {
            closest_t = ts[0];
            closest_sphere = spheres[i];
        }
        if (ts[1] < closest_t && min_t < ts[1] && ts[1] < max_t) {
            closest_t = ts[1];
            closest_sphere = spheres[i];
        }
    }

    if (closest_t == INF) {
        return BACKGROUND_COLOR;
    }

    uint8_t* color = closest_sphere.color;
    return color;
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    float camera[3] = {0, 0, 0};

    // Define scene
    Sphere spheres[3] = {
        Sphere({0, -1.0f, 3.0f}, -1.0f, {255, 0, 0}),
        Sphere({2.0f, 0, 4.0f}, 1.0f, {0, 0, 255}),
        Sphere({-2.0f, 0, 4.0f}, 1.0f, {0, 255, 0})
    };

    // Fill up the viewport with white pixels. Test that our conversions are
    // correct, from viewport to canvas.
    uint8_t color[3] = {255, 255, 255};
    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float* direction = canvas_to_viewport(x, y, width, height);
            const uint8_t* color = trace_ray(camera, direction, 1.0f, INF, spheres, 3);
            put_pixel(data, width, height, x, y, color);
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}