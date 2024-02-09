/*
Raycast 01
==========
Implements a simple raycast in a scene with 3 spheres.

```bash
g++ raytracer-01.cpp -o main.out -std=c++20
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/02-basic-raytracing.html
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;
const rgb BACKGROUND_COLOR = {255, 255, 255};

// Canvas

bool put_pixel(
    uint8_t data[][3],
    int32_t width,
    int32_t height,
    int32_t x,
    int32_t y,
    const rgb color
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

// Linear Algebra

// Compute dot product between two 3d vectors
float dot(float3 v1, float3 v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Elementwise subtraction between two 3d vectors
float3 subtract(float3 a, float3 b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

// Ray tracing

struct Sphere {
    float3 center;
    float radius;
    rgb color;

    Sphere() {}

    Sphere(const float3& v_center, float v_radius, const rgb& v_color) {
        center = v_center;
        radius = v_radius;
        color = v_color;
    }
};

// Convert 2d pixel coordinates to 3d viewport coordinates.
float3 canvas_to_viewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (float) x / width, (float) y / height, 1 };
}


// Computes intersection of ray with spheres. Returns solutions in terms of
// line parameter t.
std::array<float, 2> intersect_ray_with_sphere(
    float3 origin,
    float3 direction,
    Sphere sphere
) {
    float3 center = subtract(origin, sphere.center);

    float a = dot(direction, direction);
    float b = 2 * dot(center, direction);
    float c = dot(center, center) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant > 0) {
        return {
            (-b + sqrt(discriminant)) / (2 * a),
            (-b - sqrt(discriminant)) / (2 * a)
        };
    }
    return {INFINITY, INFINITY};
}

// Traces a ray against the spheres in the scene
rgb trace_ray(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    std::vector<Sphere> spheres
) {
    float closest_t = INFINITY;
    Sphere closest_sphere;

    for (int i = 0; i < spheres.size(); i++) {
        std::array<float, 2> ts = intersect_ray_with_sphere(origin, direction, spheres[i]);
        if (ts[0] < closest_t && min_t < ts[0] && ts[0] < max_t) {
            closest_t = ts[0];
            closest_sphere = spheres[i];
        }
        if (ts[1] < closest_t && min_t < ts[1] && ts[1] < max_t) {
            closest_t = ts[1];
            closest_sphere = spheres[i];
        }
    }

    if (closest_t == INFINITY) {
        return BACKGROUND_COLOR;
    }

    rgb color = closest_sphere.color;
    return color;
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    float3 camera = {0, 0, 0};

    // Define scene
    std::vector<Sphere> spheres = {
        Sphere({0, -1.0f, 3.0f}, -1.0f, {255, 0, 0}),
        Sphere({2.0f, 0, 4.0f}, 1.0f, {0, 0, 255}),
        Sphere({-2.0f, 0, 4.0f}, 1.0f, {0, 255, 0})
    };

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float3 direction = canvas_to_viewport(x, y, width, height);
            rgb color = trace_ray(camera, direction, 1.0f, INFINITY, spheres);
            put_pixel(data, width, height, x, y, color);
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}