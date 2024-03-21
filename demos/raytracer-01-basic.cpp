/*
Raycast 01
==========
Implements a simple raytracer in a scene with 3 spheres.

In theory, we'd need to use double precision to match Javascript's default fp64
number format, per MDN's web docs; this turns out to be necessary later on to
avoid shadow acne in Raytracer 04:
https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Number#number_encoding

However, it makes the C++ implementation here run much slower (up to 40%) than
the corresponding Javascript program.

Timing: 20ms
JS: ~120ms

```bash
g++ raytracer-01-basic.cpp -o main.out -std=c++20 -Ofast
OUT=true ./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-01.html
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;

// Canvas

bool PutPixel(
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
float DotProduct(float3 v1, float3 v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Elementwise subtraction between two 3d vectors. First minus second.
float3 Subtract(float3 v1, float3 v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
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

struct Scene {
    std::vector<Sphere> spheres;
    rgb background_color;

    Scene(std::vector<Sphere> v_spheres, rgb v_background_color) {
        spheres = v_spheres;
        background_color = v_background_color;
    }
};

// Convert 2d pixel coordinates to 3d viewport coordinates.
float3 CanvasToViewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (float) x / width, (float) y / height, 1 };
}

// Computes intersection of ray with spheres. Returns solutions in terms of
// line parameter t.
std::vector<float> IntersectRaySphere(
    float3 origin,
    float3 direction,
    Sphere sphere
) {
    float3 difference = Subtract(origin, sphere.center);

    float a = DotProduct(direction, direction);
    float b = 2 * DotProduct(difference, direction);
    float c = DotProduct(difference, difference) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {INFINITY, INFINITY};
    }

    return {
        (-b + sqrt(discriminant)) / (2 * a),
        (-b - sqrt(discriminant)) / (2 * a)
    };
}

// Traces a ray against the spheres in the scene
rgb TraceRay(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    Scene scene
) {
    float closest_t = INFINITY;
    Sphere closest_sphere;

    for (int i = 0; i < scene.spheres.size(); i++) {
        std::vector<float> ts = IntersectRaySphere(origin, direction, scene.spheres[i]);
        if (ts[0] < closest_t && min_t < ts[0] && ts[0] < max_t) {
            closest_t = ts[0];
            closest_sphere = scene.spheres[i];
        }
        if (ts[1] < closest_t && min_t < ts[1] && ts[1] < max_t) {
            closest_t = ts[1];
            closest_sphere = scene.spheres[i];
        }
    }

    if (closest_t == INFINITY) {
        return scene.background_color;
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
        Sphere({0, -1, 3}, -1, {255, 0, 0}),
        Sphere({2, 0, 4}, 1, {0, 0, 255}),
        Sphere({-2, 0, 4}, 1, {0, 255, 0})
    };
    Scene scene = Scene(spheres, {255, 255, 255});

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float3 direction = CanvasToViewport(x, y, width, height);
            rgb color = TraceRay(camera, direction, 1, INFINITY, scene);
            PutPixel(data, width, height, x, y, color);
        }
    }

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}