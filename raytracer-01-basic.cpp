/*
Raycast 01
==========
Implements a simple raycast in a scene with 3 spheres.

Note the use of double precision to match Javascript's default fp64 number
format, per MDN's web docs; this turns out to be necessary later on to avoid
shadow acne in Raytracer 04:
https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Number#number_encoding

Timing: 30ms+
JS: 115ms+

```bash
g++ raytracer-01-basic.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-01.html
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<double, 3> double3;
typedef std::array<uint8_t, 3> rgb;

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
double dot_product(double3 v1, double3 v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Elementwise subtraction between two 3d vectors. First minus second.
double3 subtract(double3 v1, double3 v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

// Ray tracing

struct Sphere {
    double3 center;
    double radius;
    rgb color;

    Sphere() {}

    Sphere(const double3& v_center, double v_radius, const rgb& v_color) {
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
double3 canvas_to_viewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (double) x / width, (double) y / height, 1 };
}

// Computes intersection of ray with spheres. Returns solutions in terms of
// line parameter t.
std::vector<double> intersect_ray_with_sphere(
    double3 origin,
    double3 direction,
    Sphere sphere
) {
    double3 difference = subtract(origin, sphere.center);

    double a = dot_product(direction, direction);
    double b = 2 * dot_product(difference, direction);
    double c = dot_product(difference, difference) - sphere.radius * sphere.radius;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {INFINITY, INFINITY};
    }

    return {
        (-b + sqrt(discriminant)) / (2 * a),
        (-b - sqrt(discriminant)) / (2 * a)
    };
}

// Traces a ray against the spheres in the scene
rgb trace_ray(
    double3 origin,
    double3 direction,
    double min_t,
    double max_t,
    Scene scene
) {
    double closest_t = INFINITY;
    Sphere closest_sphere;

    for (int i = 0; i < scene.spheres.size(); i++) {
        std::vector<double> ts = intersect_ray_with_sphere(origin, direction, scene.spheres[i]);
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
    double3 camera = {0, 0, 0};

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
            double3 direction = canvas_to_viewport(x, y, width, height);
            rgb color = trace_ray(camera, direction, 1, INFINITY, scene);
            put_pixel(data, width, height, x, y, color);
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}