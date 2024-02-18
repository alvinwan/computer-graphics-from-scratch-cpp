/*
Raycast 05
==========
Adds reflections

Timing: 269ms
JS: ~280ms

```bash
g++ raytracer-05-reflection.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-05.html
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;
const float EPSILON = 0.001;

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

// Length of vector
float Length(float3 vec) {
    return sqrt(DotProduct(vec, vec));
}

// Broadcasted Multiply between scalar and a vector
float3 Multiply(float k, float3 vec) {
    return {k * vec[0], k * vec[1], k * vec[2]};
}

// Elementwise addition between two 3d vectors
float3 Add(float3 v1, float3 v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

// Elementwise subtraction between two 3d vectors. First minus second.
float3 Subtract(float3 v1, float3 v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

rgb Clamp(float3 vec) {
    return {
        (uint8_t) std::round(std::clamp<float>(vec[0], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(vec[1], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(vec[2], 0, 255)),
    };
}

// Ray tracing

struct Sphere {
    float3 center;
    float radius;
    float3 color;
    float specular;
    float reflective;

    Sphere() {}

    Sphere(const float3& v_center, float v_radius, const float3& v_color, float v_specular, float v_reflective) {
        center = v_center;
        radius = v_radius;
        color = v_color;
        specular = v_specular;
        reflective = v_reflective;
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
    float3 background_color;

    Scene(std::vector<Sphere> v_spheres, std::vector<Light> v_lights, float3 v_background_color) {
        spheres = v_spheres;
        lights = v_lights;
        background_color = v_background_color;
    }
};

// Convert 2d pixel coordinates to 3d viewport coordinates.
float3 CanvasToViewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (float) x / width, (float) y / height, 1 };
}

// Computes intersection of ray with spheres. Returns solutions in terms of
// line parameter t.
std::vector<float> intersect_ray_with_sphere(
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

// Find the closest intersection between a ray and the spheres in the scene.
std::tuple<Sphere, float> ClosestIntersection(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    Scene scene
) {
    float closest_t = INFINITY;
    Sphere closest_sphere;

    for (int i = 0; i < scene.spheres.size(); i++) {
        std::vector<float> ts = intersect_ray_with_sphere(origin, direction, scene.spheres[i]);
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

// Compute the reflection of a ray on a surface defined by its normal
float3 ReflectRay(float3 ray, float3 normal) {
    return Subtract(Multiply(2 * DotProduct(ray, normal), normal), ray);
}

// Compute lighting for the scene
float ComputeLighting(float3 point, float3 normal, float3 view, float specular, Scene scene) {
    float intensity = 0;
    if (abs(Length(normal) - 1) > 0.0001) {
        std::cerr << "Error: Normal is not length 1 (" << Length(normal) << ")" << std::endl;
        return INFINITY;
    }

    float length_v = Length(view);

    for (int i = 0; i < scene.lights.size(); i++) {
        Light light = scene.lights[i];
        if (light.ltype == AMBIENT) {
            intensity += light.intensity;
        } else {
            float3 vec_l;
            float shadow_t_max;

            if (light.ltype == POINT) {
                vec_l = Subtract(light.position, point);
                shadow_t_max = 1;
            } else {  // Light.DIRECTIONAL
                vec_l = light.position;
                shadow_t_max = INFINITY;
            }

            // Shadow check
            std::tuple<Sphere, float> intersection = ClosestIntersection(point, vec_l, EPSILON, shadow_t_max, scene);
            if (std::get<1>(intersection) != INFINITY) continue;

            // Diffuse
            float n_dot_l = DotProduct(normal, vec_l);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (Length(vec_l));
            }

            // Specular, where vec_r is the 'perfect' reflection ray
            if (specular != -1) {
                float3 vec_r = ReflectRay(vec_l, normal);
                float r_dot_v = DotProduct(vec_r, view);
                if (r_dot_v > 0) {
                    intensity += light.intensity * pow(r_dot_v / (Length(vec_r) * length_v), specular);
                }
            }
        }
    }

    return intensity;
}

// Traces a ray against the spheres in the scene
float3 TraceRay(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    int32_t recursion_depth,
    Scene scene
) {
    std::tuple<Sphere, float> intersection = ClosestIntersection(origin, direction, min_t, max_t, scene);
    Sphere closest_sphere = std::get<0>(intersection);
    float closest_t = std::get<1>(intersection);

    if (closest_t == INFINITY) {
        return scene.background_color;
    }

    float3 point = Add(origin, Multiply(closest_t, direction));
    float3 normal = Subtract(point, closest_sphere.center);
    normal = Multiply(1 / Length(normal), normal);

    float intensity = ComputeLighting(point, normal, Multiply(-1, direction), closest_sphere.specular, scene);
    float3 local_color = Multiply(intensity, closest_sphere.color);

    // If we hit the recursion limit or the sphere is not reflective, finish
    float reflective = closest_sphere.reflective;
    if (recursion_depth <= 0 || reflective <= 0) {
        return local_color;
    }

    // Compute the reflected color
    float3 reflected_ray = ReflectRay(Multiply(-1, direction), normal);
    // NOTE: Below addresses 'shadow' acne (i.e., flickering shadow), which you
    // can see here: https://imgur.com/a/ycB69zX. To fix this, move the starting
    // point of the reflection ray along the normal, to prevent self-
    // intersection. JS demos don't have this problem because JS uses FP64 by
    // default. Additionally, other spheres don't have this problem because only
    // the biggest radius=5000 one has reflection rays that are near parallel.
    point = Add(point, Multiply(EPSILON, normal));
    const float3 reflected_color = TraceRay(point, reflected_ray, EPSILON, INFINITY, recursion_depth - 1, scene);

    return Add(Multiply(1 - reflective, local_color), Multiply(reflective, reflected_color));
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    float3 camera = {0, 0, 0};

    // Define scene
    std::vector<Sphere> spheres = {
        Sphere({0, -1, 3}, 1, {255, 0, 0}, 500, 0.2),
        Sphere({-2, 0, 4}, 1, {0, 255, 0}, 10, 0.4),
        Sphere({2, 0, 4}, 1, {0, 0, 255}, 500, 0.3),
        Sphere({0, -5001, 0}, 5000, {255, 255, 0}, 1000, 0.5)
    };
    std::vector<Light> lights = {
        Light(AMBIENT, 0.2, {0, 0, 0}),
        Light(POINT, 0.6, {2, 1, 0}),
        Light(DIRECTIONAL, 0.2, {1, 4, 4})
    };
    Scene scene = Scene(spheres, lights, {0, 0, 0});

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float3 direction = CanvasToViewport(x, y, width, height);
            float3 color = TraceRay(camera, direction, 1, INFINITY, 3, scene);
            PutPixel(data, width, height, x, y, Clamp(color));
        }
    }

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}