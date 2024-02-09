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

const float INF = std::numeric_limits<float>::infinity();
const uint8_t BACKGROUND_COLOR[3] = {255, 255, 255};

// Canvas

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

// Linear Algebra

// Compute dot product between two 3d vectors
float dot(float v1[3], float v2[3]) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Length of vector
float length(float vec[3]) {
    return sqrt(dot(vec, vec));
}

// Broadcasted multiply between vector and scalar
float* multiply(float k, float vec[3]) {
    float* point = new float[3];
    point[0] = k * vec[0];
    point[1] = k * vec[1];
    point[2] = k * vec[2];
    return point;
}

// Elementwise addition between two 3d vectors
float* add(float v1[3], float v2[3]) {
    float* point = new float[3];
    point[0] = v1[0] + v2[0];
    point[1] = v1[1] + v2[1];
    point[2] = v1[2] + v2[2];
    return point;
}

// Elementwise subtraction between two 3d vectors
float* subtract(float v1[3], float v2[3]) {
    float* point = new float[3];
    point[0] = v1[0] - v2[0];
    point[1] = v1[1] - v2[1];
    point[2] = v1[2] - v2[2];
    return point;
}

uint8_t* raw_to_color(float vec[3]) {
    uint8_t* color = new uint8_t[3];
    color[0] = (uint8_t) std::round(std::clamp<float>(vec[0], 0.0f, 255.0f));
    color[1] = (uint8_t) std::round(std::clamp<float>(vec[1], 0.0f, 255.0f));
    color[2] = (uint8_t) std::round(std::clamp<float>(vec[2], 0.0f, 255.0f));
    return color;
}

// Ray tracing

struct Sphere {
    float center[3];
    float radius;
    float color[3];
    float specular;

    Sphere() {}

    Sphere(const std::initializer_list<float>& c, float r, const std::initializer_list<uint8_t>& col, float s) {
        std::copy(c.begin(), c.end(), center);
        radius = r;
        std::copy(col.begin(), col.end(), color);
        specular = s;
    }
};

enum LightType {AMBIENT, POINT, DIRECTIONAL};

struct Light {
    LightType ltype;
    float intensity;
    float position[3];

    Light() {}

    Light(LightType l, float i, const std::initializer_list<float>& p) {
        ltype = l;
        intensity = i;
        std::copy(p.begin(), p.end(), position);
    }
};

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
    if (discriminant >= 0) {
        out[0] = (-b + sqrt(discriminant)) / (2 * a);
        out[1] = (-b - sqrt(discriminant)) / (2 * a);
    } else {
        out[0] = INF;
        out[1] = INF;
    }
    return out;
}

// Find the closest intersection between a ray and the spheres in the scene.
std::tuple<Sphere, float> closest_intersection(
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

    return std::make_tuple(closest_sphere, closest_t);
}

// Compute lighting for the scene
float compute_lighting(
    float point[3],
    float normal[3],
    float view[3],
    float specular,
    Sphere spheres[],
    int32_t num_spheres,
    Light lights[],
    int32_t num_lights
) {
    float intensity = 0;
    if (abs(length(normal) - 1.0f) > 0.0001f) {
        std::cerr << "Error: Normal is not length 1 (" << length(normal) << ")" << std::endl;
        return INF;
    }

    float length_v = length(view);

    for (int i = 0; i < num_lights; i++) {
        Light light = lights[i];
        if (light.ltype == AMBIENT) {
            intensity += light.intensity;
        } else {
            float* vec_l;
            float shadow_t_max;
            if (light.ltype == POINT) {
                vec_l = subtract(light.position, point);
                shadow_t_max = 1.0f;
            } else {  // Light.DIRECTIONAL
                vec_l = light.position;
                shadow_t_max = INF;
            }

            // Shadow check
            std::tuple<Sphere, float> intersection = closest_intersection(point, vec_l, 0.0001, shadow_t_max, spheres, num_spheres);
            if (std::get<1>(intersection) != INF) continue;

            // Diffuse
            float n_dot_l = dot(normal, vec_l);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (length(vec_l));
            }

            // Specular, where vec_r is the 'perfect' reflection ray
            if (specular != -1) {
                float* vec_r = subtract(multiply(2 * dot(normal, vec_l), normal), vec_l);
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
const uint8_t* trace_ray(
    float origin[3],
    float direction[3],
    float min_t,
    float max_t,
    Sphere spheres[],
    int32_t num_spheres,
    Light lights[],
    int32_t num_lights
) {
    std::tuple<Sphere, float> intersection = closest_intersection(origin, direction, min_t, max_t, spheres, num_spheres);
    Sphere closest_sphere = std::get<0>(intersection);
    float closest_t = std::get<1>(intersection);

    if (closest_t == INF) {
        return BACKGROUND_COLOR;
    }

    float* point = add(origin, multiply(closest_t, direction));
    float* normal = subtract(point, closest_sphere.center);
    normal = multiply(1.0f / length(normal), normal);

    float intensity = compute_lighting(point, normal, multiply(-1, direction), closest_sphere.specular, spheres, num_spheres, lights, num_lights);
    float* raw = multiply(intensity, closest_sphere.color);
    uint8_t* color = raw_to_color(raw);

    return color;
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    float camera[3] = {0, 0, 0};

    // Define scene
    Sphere spheres[4] = {
        Sphere({0, -1.0f, 3.0f}, -1.0f, {255, 0, 0}, 500.0f),
        Sphere({2.0f, 0, 4.0f}, 1.0f, {0, 0, 255}, 500.0f),
        Sphere({-2.0f, 0, 4.0f}, 1.0f, {0, 255, 0}, 10.0f),
        Sphere({0, -5001.0f, 0}, 5000.0f, {255, 255, 0}, 1000.0f)
    };
    Light lights[3] = {
        Light(AMBIENT, 0.2f, {0, 0, 0}),
        Light(POINT, 0.6f, {2, 1, 0}),
        Light(DIRECTIONAL, 0.2f, {1, 4, 4})
    };

    uint8_t color[3] = {255, 255, 255};
    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float* direction = canvas_to_viewport(x, y, width, height);
            const uint8_t* color = trace_ray(camera, direction, 1.0f, INF, spheres, 4, lights, 3);
            put_pixel(data, width, height, x, y, color);
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}