/*
Raycast 08 - Constructive Solid Geometry
========================================
Adds support for constructive solid geometry

Timing: 1.32s

```bash
g++ raytracer-08-csg.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

The original "Computer Graphics from Scratch" book doesn't contain a demo for
this additional functionality. The original description of the algorithm is in
the "Extending the Raytracer" chapter:
https://gabrielgambetta.com/computer-graphics-from-scratch/05-extending-the-raytracer.html#constructive-solid-geometry
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<float, 3> float3;
typedef std::array<float3, 3> float33;
typedef std::array<uint8_t, 3> rgb;
const float EPSILON = 0.001;

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
float dot_product(float3 v1, float3 v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Length of vector
float length(float3 vec) {
    return sqrt(dot_product(vec, vec));
}

// Broadcasted multiply between scalar and a vector
float3 multiply(float k, float3 vec) {
    return {k * vec[0], k * vec[1], k * vec[2]};
}

// Elementwise addition between two 3d vectors
float3 add(float3 v1, float3 v2) {
    return {v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]};
}

// Elementwise subtraction between two 3d vectors. First minus second.
float3 subtract(float3 v1, float3 v2) {
    return {v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]};
}

// Multiply 3x3 matrix and 3x1 vector
float3 matmul(float33 matrix, float3 vector) {
    float3 out = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i] += vector[j] * matrix[i][j];
        }
    }
    return out;
}

// Implements cross product
float3 cross_product(float3 v1, float3 v2) {
    return {
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0]
    };
}

rgb clamp(float3 vec) {
    return {
        (uint8_t) std::round(std::clamp<float>(vec[0], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(vec[1], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(vec[2], 0, 255)),
    };
}

// Ray tracing

struct Material {
    float3 color;
    float specular;
    float reflective;

    Material() {}

    Material(const float3& v_color, float v_specular, float v_reflective) {
        color = v_color;
        specular = v_specular;
        reflective = v_reflective;
    }
};

struct Object {
    Material material;

    Object() {}

    Object(Material v_material) {
        material = v_material;
    }
};

struct Sphere : Object {
    float3 center;
    float radius;

    Sphere() {}

    Sphere(const float3& v_center, float v_radius, Material v_material)
    : Object(v_material) {
        center = v_center;
        radius = v_radius;
    }
};

// Computes intersection of a ray with sphere. Returns solution in terms of
// line parameter t.
std::vector<float> intersect(Sphere sphere, float3 origin, float3 direction) {
    float3 difference = subtract(origin, sphere.center);

    float a = dot_product(direction, direction);
    float b = 2 * dot_product(difference, direction);
    float c = dot_product(difference, difference) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {INFINITY, INFINITY};
    }

    float t0 = (-b + sqrt(discriminant)) / (2 * a);
    float t1 = (-b - sqrt(discriminant)) / (2 * a);
    return {std::min(t0, t1), std::max(t0, t1)};
}

// Returns normal for a given point on the surface of a sphere.
float3 get_normal_of(Sphere sphere, float3 point) {
    float3 normal = subtract(point, sphere.center);
    normal = multiply(1 / length(normal), normal);
    return normal;
}

// Checks if the point lies on this sphere's surface.
bool has_point(Sphere sphere, float3 point) {
    return std::abs(length(subtract(point, sphere.center)) - sphere.radius) < EPSILON;
}

struct Plane : Object {
    float3 normal;
    float distance;

    Plane() {}

    Plane(float3 v_normal, float v_distance, Material v_material)
    : Object(v_material) {
        normal = v_normal;
        distance = v_distance;
    }
};

// Computes intersection of a ray with plane. Returns solution in terms of
// line parameter t.
std::vector<float> intersect(Plane plane, float3 origin, float3 direction) {
    float denominator = dot_product(plane.normal, direction);
    // Plane is parallel to ray (=0) or we're looking at back side (<0)
    if (denominator >= 0) return {INFINITY};

    float t = -(plane.distance + dot_product(plane.normal, origin)) / denominator;
    if (t < 0) return {INFINITY}; // Triangle is 'behind' the ray
    return {t};
}

struct Triangle : Object {
    float3 a;
    float3 b;
    float3 c;
    Plane plane;

    Triangle(float3 v_a, float3 v_b, float3 v_c, Material v_material)
    : Object(v_material) {
        a = v_a;
        b = v_b;
        c = v_c;

        // Computes plane containing triangle. Note this is redundant
        // with vertex information but allows us to cache computation.
        float3 ab = subtract(b, a);
        float3 ac = subtract(c, a);
        float3 normal = cross_product(ab, ac);
        float magnitude = length(normal);
        if (magnitude == 0) {
            std::cerr << "Error: Triangle is degenerate (i.e., points are collinear)" << std::endl;
            return;
        }
        normal = multiply(1 / magnitude, normal);

        float distance = -dot_product(normal, a);
        plane = Plane(normal, distance, v_material);
    }
};

// Compute on which side a point lies, relative to the line defined by two
// points. This is the sign that you would get by computing
// normal.dot_product(point) - distance.
float sign(float3 point, float3 a, float3 b, float3 normal) {
    float3 candidate = subtract(point, a);
    float3 edge = subtract(b, a);
    return dot_product(normal, cross_product(candidate, edge));
}

// Computes intersection of a ray with triangle. Returns solution in terms
// of line parameter t.
std::vector<float> intersect(Triangle triangle, float3 origin, float3 direction) {
    // Find intersection between ray and plane
    float t = intersect(triangle.plane, origin, direction)[0];
    if (t == INFINITY) return {INFINITY};
    float3 point = add(origin, multiply(t, direction));

    // Check if the intersection lies in the triangle. NOTE the ordering of
    // points must be consistently clockwise OR counter clockwise, in this
    // calculation. If *any sign is negative, ray does not intersect the
    // triangle.
    if (sign(point, triangle.a, triangle.b, triangle.plane.normal) >= 0) return {INFINITY};
    if (sign(point, triangle.b, triangle.c, triangle.plane.normal) >= 0) return {INFINITY};
    if (sign(point, triangle.c, triangle.a, triangle.plane.normal) >= 0) return {INFINITY};

    // The intersection point lies in the 'correct' half plane for every
    // edge, so it lies in the triangle.
    return {t};
}

enum Operation { AND, OR, MINUS };

struct CSG : Object {
    Sphere sphere1;
    Sphere sphere2;
    Operation operation;

    CSG() {}

    CSG(Sphere v_sphere1, Sphere v_sphere2, Operation v_operation, Material v_material)
    : Object(v_material) {
        sphere1 = v_sphere1;
        sphere2 = v_sphere2;
        operation = v_operation;
    }
};

std::vector<float> intersect(CSG csg, float3 origin, float3 direction) {
    std::vector<float> ts1 = intersect(csg.sphere1, origin, direction);
    std::vector<float> ts2 = intersect(csg.sphere2, origin, direction);

    // Objects without volume will only return a single intersection, such
    // as a plane or triangle. Objects *with volume will return an even
    // number of interesections.
    if (ts1.size() % 2 != 0 || ts2.size() % 2 != 0) {
        std::cerr << "Error: Both objects in a CSG must have volume" << std::endl;
        return {INFINITY, INFINITY};
    }

    // Each operation has a different response, depending on which of the
    // two objects the ray intersects.
    if (csg.operation == AND) {
        // In an AND operation, *both volumes must be intersected
        if (ts1[0] != INFINITY && ts2[0] != INFINITY) {
            return {std::min(ts1[1], ts2[1]), std::max(ts1[0], ts2[0])};
        }
        return {INFINITY, INFINITY};
    } else if (csg.operation == OR) {
        if (ts1[0] != INFINITY && ts2[0] != INFINITY) {
            return {std::min(ts1[0], ts2[0]), std::max(ts1[1], ts2[1])};
        }
        // In an OR operation, AT LEAST one volume needs to be intersected
        if (ts1[0] == INFINITY || ts2[0] == INFINITY) {
            return {std::min(ts1[0], ts2[0]), std::min(ts1[1], ts2[1])};
        }
        return {INFINITY, INFINITY}; // Intersects neither
    } else {  // operation == MINUS
        // In a MINUS operation, the first object must be intersected.
        // If you order all the t's together, we allow BABA, ABAB, AA, ABBA,
        // and disallow BB, BAAB. We can check this by ensuring that A[0]
        // is in front OR A[1] is in the back.
        if (ts1[0] != INFINITY && (ts1[0] < ts2[0] || ts1[1] > ts2[1])) {
            return {std::min(ts1[0], ts2[0]), std::max(ts1[1], ts2[1])};
        }
        return {INFINITY, INFINITY};
    }
}

float3 get_normal_of(CSG csg, float3 point) {
    if (has_point(csg.sphere1, point)) {
        // No matter the operation, if the point of intersection lies on the
        // first object, always use the first object's normal, normally.
        return get_normal_of(csg.sphere1, point);
    }
    // Otherwise, the point lies on the second object. In this case, we need
    // to handle operation-by-operation.
    float3 normal = get_normal_of(csg.sphere2, point);
    if (csg.operation == AND || csg.operation == OR) {
        // If operation is AND or OR, use the normal, normally.
        return normal;
    }
    // If the operation is a MINUS, and the point of intersection lies on
    // the second object, the normal needs to be flipped so it's on the
    // inside of object2.
    return multiply(-1, normal);
}

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

struct Camera {
    float3 position;
    float33 rotation;

    Camera(float3 v_position, float33 v_rotation) {
        position = v_position;
        rotation = v_rotation;
    }
};

struct Scene {
    std::vector<Sphere> spheres;
    std::vector<Triangle> triangles;
    std::vector<CSG> csgs;
    std::vector<Light> lights;
    float3 background_color;

    Scene(std::vector<Sphere> v_spheres, std::vector<Triangle> v_triangles, std::vector<CSG> v_csgs, std::vector<Light> v_lights, float3 v_background_color) {
        spheres = v_spheres;
        triangles = v_triangles;
        csgs = v_csgs;
        lights = v_lights;
        background_color = v_background_color;
    }
};

// Convert 2d pixel coordinates to 3d viewport coordinates.
float3 canvas_to_viewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (float) x / width, (float) y / height, 1 };
}

// Holds intersection information for a raycast
struct Intersection {
    Material material;
    float3 point;
    float3 normal;
    bool is_valid;

    Intersection() {
        is_valid = false;
    }

    Intersection(Material v_material, float3 v_point, float3 v_normal) {
        material = v_material;
        point = v_point;
        normal = v_normal;
        is_valid = true;
    }
};

// Find the closest intersection between a ray and the spheres in the scene.
Intersection closest_intersection(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    Scene scene
) {
    float closest_t = INFINITY;
    Material material;
    float3 normal;

    for (int i = 0; i < scene.spheres.size(); i++) {
        Sphere sphere = scene.spheres[i];

        std::vector<float> ts = intersect(sphere, origin, direction);
        for (int j = 0; j < ts.size(); j++) {
            if (ts[j] < closest_t && min_t < ts[j] && ts[j] < max_t) {
                closest_t = ts[j];
                material = sphere.material;

                float3 point = add(origin, multiply(closest_t, direction));
                normal = get_normal_of(sphere, point);
            }
        }
    }

    // NOTE: Duplicated code for triangles, to use compile-time polymorphism
    for (int i = 0; i < scene.triangles.size(); i++) {
        Triangle triangle = scene.triangles[i];

        float t = intersect(triangle, origin, direction)[0];
        if (t < closest_t && min_t < t && t < max_t) {
            closest_t = t;
            material = triangle.material;
            normal = triangle.plane.normal;
        }
    }

    // NOTE: Duplicated code for CSGs, to use compile-time polymorphism
    for (int i = 0; i < scene.csgs.size(); i++) {
        CSG csg = scene.csgs[i];

        std::vector<float> ts = intersect(csg, origin, direction);
        for (int j = 0; j < ts.size(); j++) {
            if (ts[j] < closest_t && min_t < ts[j] && ts[j] < max_t) {
                closest_t = ts[j];
                material = csg.material;

                float3 point = add(origin, multiply(closest_t, direction));
                normal = get_normal_of(csg, point);
            }
        }
    }

    // sets intersection.is_valid=false
    if (closest_t == INFINITY) return Intersection();

    // sets intersection.is_valid = true
    float3 point = add(origin, multiply(closest_t, direction));
    return Intersection(material, point, normal);
}

// Compute the reflection of a ray on a surface defined by its normal
float3 reflect_ray(float3 ray, float3 normal) {
    return subtract(multiply(2 * dot_product(ray, normal), normal), ray);
}

// Compute lighting for the scene
float compute_lighting(float3 point, float3 normal, float3 view, float specular, Scene scene) {
    float intensity = 0;
    if (abs(length(normal) - 1) > 0.0001) {
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
                shadow_t_max = 1;
            } else {  // Light.DIRECTIONAL
                vec_l = light.position;
                shadow_t_max = INFINITY;
            }

            // Shadow check
            Intersection intersection = closest_intersection(point, vec_l, EPSILON, shadow_t_max, scene);
            if (intersection.is_valid) continue;

            // Diffuse
            float n_dot_l = dot_product(normal, vec_l);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (length(vec_l));
            }

            // Specular, where vec_r is the 'perfect' reflection ray
            if (specular != -1) {
                float3 vec_r = reflect_ray(vec_l, normal);
                float r_dot_v = dot_product(vec_r, view);
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
    int32_t recursion_depth,
    Scene scene
) {
    Intersection intersection = closest_intersection(origin, direction, min_t, max_t, scene);
    if (!intersection.is_valid) {
        return scene.background_color;
    }

    float intensity = compute_lighting(
        intersection.point,
        intersection.normal,
        multiply(-1, direction),
        intersection.material.specular,
        scene);
    float3 local_color = multiply(intensity, intersection.material.color);

    // If we hit the recursion limit or the sphere is not reflective, finish
    float reflective = intersection.material.reflective;
    if (recursion_depth <= 0 || reflective <= 0) {
        return local_color;
    }

    // Compute the reflected color
    float3 reflected_ray = reflect_ray(multiply(-1, direction), intersection.normal);
    const float3 reflected_color = trace_ray(intersection.point, reflected_ray, EPSILON, INFINITY, recursion_depth - 1, scene);

    return add(multiply(1 - reflective, local_color), multiply(reflective, reflected_color));
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    float3 position = {3, 0, 1};
    float33 rotation = {{
        {{0.7071, 0, -0.7071}},
        {{     0, 1,       0}},
        {{0.7071, 0,  0.7071}}
    }};
    Camera camera = Camera(position, rotation);

    // Define scene
    std::vector<Sphere> spheres = {
        Sphere({0, -5001, 0}, 5000, Material({255, 255, 0}, 1000, 0.5)),
    };
    std::vector<Triangle> triangles = {
        Triangle({2, 0, 6}, {-2, 0, 6}, {0, 2, 4}, Material({0, 255, 255}, 500, 0.4))
    };
    std::vector<CSG> csgs = {
        CSG(
            Sphere({-2, 0, 4}, 1, Material()),
            Sphere({-2, 1, 4}, 1, Material()),
            OR, Material({0, 255, 0}, 10, 0.4)
        ),
        CSG(
            Sphere({0, -1, 3}, 1, Material()),
            Sphere({0, 0, 3}, 1, Material()),
            AND, Material({255, 0, 0}, 500, 0.2f)
        ),
        CSG(
            Sphere({2, 0, 4}, 1, Material()),
            Sphere({2, 1, 3}, 1, Material()),
            MINUS, Material({0, 0, 255}, 500, 0.3)
        ),
    };

    std::vector<Light> lights = {
        Light(AMBIENT, 0.2, {0, 0, 0}),
        Light(POINT, 0.6, {2, 1, 0}),
        Light(DIRECTIONAL, 0.2, {1, 4, 4})
    };
    Scene scene = Scene(spheres, triangles, csgs, lights, {0, 0, 0});

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float3 direction = matmul(camera.rotation, canvas_to_viewport(x, y, width, height));
            float3 color = trace_ray(camera.position, direction, 1, INFINITY, 3, scene);
            put_pixel(data, width, height, x, y, clamp(color));
        }
    }

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}