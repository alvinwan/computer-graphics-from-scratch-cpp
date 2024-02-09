/*
Raycast 04
==========
Adds shadows

```bash
g++ raytracer-04.cpp -o main.out -std=c++20
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-04.html
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;
const float3 BACKGROUND_COLOR = {255, 255, 255};
const float EPSILON = 0.0001f;

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

// Length of vector
float length(float3 vec) {
    return sqrt(dot(vec, vec));
}

// Broadcasted multiply between scalar and a vector
float3 multiply(float k, float3 vec) {
    return {k * vec[0], k * vec[1], k * vec[2]};
}

// Elementwise addition between two 3d vectors
float3 add(float3 a, float3 b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

// Elementwise subtraction between two 3d vectors
float3 subtract(float3 a, float3 b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

rgb clamp(float3 vec) {
    return {
        (uint8_t) std::round(std::clamp<float>(vec[0], 0.0f, 255.0)),
        (uint8_t) std::round(std::clamp<float>(vec[1], 0.0f, 255.0)),
        (uint8_t) std::round(std::clamp<float>(vec[2], 0.0f, 255.0)),
    };
}

// Ray tracing

struct Sphere {
    float3 center;
    float radius;
    float3 color;
    float specular;

    Sphere() {}

    Sphere(const float3& v_center, float v_radius, const float3& v_color, float v_specular) {
        center = v_center;
        radius = v_radius;
        color = v_color;
        specular = v_specular;
    }
};

enum LightType {AMBIENT, POINT, DIRECTIONAL};

struct Light {
    LightType ltype;
    float intensity;
    float3 position;

    Light() {}

    Light(LightType v_ltype, float v_intensity, const float3& v_position) {
        ltype = v_ltype;
        intensity = v_intensity;
        position = v_position;
    }
};

struct Scene {
    std::vector<Sphere> spheres;
    std::vector<Light> lights;

    Scene(std::vector<Sphere> v_spheres, std::vector<Light> v_lights) {
        spheres = v_spheres;
        lights = v_lights;
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

// Find the closest intersection between a ray and the spheres in the scene.
std::tuple<Sphere, float> closest_intersection(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    Scene scene
) {
    float closest_t = INFINITY;
    Sphere closest_sphere;

    for (int i = 0; i < scene.spheres.size(); i++) {
        std::array<float, 2> ts = intersect_ray_with_sphere(origin, direction, scene.spheres[i]);
        if (ts[0] < closest_t && min_t < ts[0] && ts[0] < max_t) {
            closest_t = ts[0];
            closest_sphere = scene.spheres[i];
        }
        if (ts[1] < closest_t && min_t < ts[1] && ts[1] < max_t) {
            closest_t = ts[1];
            closest_sphere = scene.spheres[i];
        }
    }

    return std::make_tuple(closest_sphere, closest_t);
}

// Compute lighting for the scene
float compute_lighting(float3 point, float3 normal, float3 view, float specular, Scene scene) {
    float intensity = 0;
    if (abs(length(normal) - 1.0f) > 0.0001f) {
        std::cerr << "Error: Normal is not length 1 (" << length(normal) << ")" << std::endl;
        return INFINITY;
    }

    float length_v = length(view);

    for (int i = 0; i < scene.lights.size(); i++) {
        Light light = scene.lights[i];
        if (light.ltype == AMBIENT) {
            intensity += light.intensity;
        } else {
            float3 vec_l;
            float shadow_t_max;

            if (light.ltype == POINT) {
                vec_l = subtract(light.position, point);
                shadow_t_max = 1.0f;
            } else {  // Light.DIRECTIONAL
                vec_l = light.position;
                shadow_t_max = INFINITY;
            }

            // Shadow check
            std::tuple<Sphere, float> intersection = closest_intersection(point, vec_l, EPSILON, shadow_t_max, scene);
            if (std::get<1>(intersection) != INFINITY) continue;

            // Diffuse
            float n_dot_l = dot(normal, vec_l);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (length(vec_l));
            }

            // Specular, where vec_r is the 'perfect' reflection ray
            if (specular != -1) {
                float3 vec_r = subtract(multiply(2 * dot(normal, vec_l), normal), vec_l);
                float r_dot_v = dot(vec_r, view);
                if (r_dot_v > 0) {
                    intensity += light.intensity * pow(r_dot_v / (length(vec_r) * length_v), specular);
                }
            }
        }
    }

    return intensity;
}

// Traces a ray against the spheres in the scene
float3 trace_ray(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    Scene scene
) {
    std::tuple<Sphere, float> intersection = closest_intersection(origin, direction, min_t, max_t, scene);
    Sphere closest_sphere = std::get<0>(intersection);
    float closest_t = std::get<1>(intersection);

    if (closest_t == INFINITY) {
        return BACKGROUND_COLOR;
    }

    float3 point = add(origin, multiply(closest_t, direction));
    float3 normal = subtract(point, closest_sphere.center);
    normal = multiply(1.0f / length(normal), normal);

    float intensity = compute_lighting(point, normal, multiply(-1, direction), closest_sphere.specular, scene);
    return multiply(intensity, closest_sphere.color);
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    float3 camera = {0, 0, 0};

    // Define scene
    std::vector<Sphere> spheres = {
        Sphere({0, -1.0f, 3.0f}, -1.0f, {255, 0, 0}, 500.0f),
        Sphere({2.0f, 0, 4.0f}, 1.0f, {0, 0, 255}, 500.0f),
        Sphere({-2.0f, 0, 4.0f}, 1.0f, {0, 255, 0}, 10.0f),
        Sphere({0, -5001.0f, 0}, 5000.0f, {255, 255, 0}, 1000.0f)
    };
    std::vector<Light> lights = {
        Light(AMBIENT, 0.2f, {0, 0, 0}),
        Light(POINT, 0.6f, {2, 1, 0}),
        Light(DIRECTIONAL, 0.2f, {1, 4, 4})
    };
    Scene scene = Scene(spheres, lights);

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float3 direction = canvas_to_viewport(x, y, width, height);
            float3 color = trace_ray(camera, direction, 1.0f, INFINITY, scene);
            put_pixel(data, width, height, x, y, clamp(color));
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}