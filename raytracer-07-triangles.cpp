/*
Raycast 07 - Triangles
======================
Adds support for triangle primitives to the last demo that "Computer Graphics
from Scratch" has an implementation for (Raytracer 06).

Timing: ~1.36s

```bash
g++ raytracer-07-triangles.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

The original "Computer Graphics from Scratch" book doesn't contain a demo for
this additional functionality. The original description of the algorithm is in
the "Extending the Raytracer" chapter:
https://gabrielgambetta.com/computer-graphics-from-scratch/05-extending-the-raytracer.html#other-primitives
*/
#include "bmp.h"
#include <math.h>
#include <array>

typedef std::array<double, 3> double3;
typedef std::array<double3, 3> double33;
typedef std::array<uint8_t, 3> rgb;
const double3 BACKGROUND_COLOR = {0, 0, 0};
const double EPSILON = 0.001f;

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
double dot(double3 v1, double3 v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Length of vector
double length(double3 vec) {
    return sqrt(dot(vec, vec));
}

// Broadcasted multiply between scalar and a vector
double3 multiply(double k, double3 vec) {
    return {k * vec[0], k * vec[1], k * vec[2]};
}

// Elementwise addition between two 3d vectors
double3 add(double3 a, double3 b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

// Elementwise subtraction between two 3d vectors
double3 subtract(double3 a, double3 b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

// Multiply 3x3 matrix and 3x1 vector
double3 matmul(double33 matrix, double3 vector) {
    double3 out = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i] += vector[j] * matrix[i][j];
        }
    }
    return out;
}

// Implements cross product
double3 cross_product(double3 a, double3 b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

rgb clamp(double3 vec) {
    return {
        (uint8_t) std::round(std::clamp<double>(vec[0], 0.0f, 255.0)),
        (uint8_t) std::round(std::clamp<double>(vec[1], 0.0f, 255.0)),
        (uint8_t) std::round(std::clamp<double>(vec[2], 0.0f, 255.0)),
    };
}

// Ray tracing

enum Shape { SPHERE, TRIANGLE, PLANE };

struct Object {
    Shape shape;
    double3 color;
    double specular;
    double reflective;

    Object() {}

    Object(const double3& v_color, double v_specular, double v_reflective, Shape v_shape) {
        color = v_color;
        specular = v_specular;
        reflective = v_reflective;
        shape = v_shape;
    }

    // Computes intersection of a ray with object. Returns solution in terms of
    // line parameter t.
    virtual std::vector<double> intersect(double3 origin, double3 direction) {
        return {INFINITY};
    }

    // Returns normal for a given point on the surface of this object.
    virtual double3 get_normal_of(double3 point) {
        return {0, 0, 0};
    }
};

struct Sphere : Object {
    double3 center;
    double radius;

    Sphere() {}

    Sphere(const double3& v_center, double v_radius, const double3& v_color, double v_specular, double v_reflective)
    : Object(v_color, v_specular, v_reflective, SPHERE) {
        center = v_center;
        radius = v_radius;
    }

    // Computes intersection of a ray with sphere. Returns solution in terms of
    // line parameter t.
    std::vector<double> intersect(double3 origin, double3 direction) {
        double3 offset = subtract(origin, center);

        double a = dot(direction, direction);
        double b = 2 * dot(offset, direction);
        double c = dot(offset, offset) - radius * radius;

        double discriminant = b * b - 4 * a * c;
        if (discriminant >= 0) {
            return {
                (-b + sqrt(discriminant)) / (2 * a),
                (-b - sqrt(discriminant)) / (2 * a)
            };
        }
        return {INFINITY, INFINITY};
    }

    double3 get_normal_of(double3 point) {
        double3 normal = subtract(point, center);
        normal = multiply(1.0f / length(normal), normal);
        return normal;
    }
};

struct Plane : Object {
    double3 normal;
    double distance;

    Plane() {}

    Plane(double3 v_normal, double v_distance, const double3& v_color, double v_specular, double v_reflective)
    : Object(v_color, v_specular, v_reflective, PLANE) {
        normal = v_normal;
        distance = v_distance;
    }

    // Computes intersection of a ray with plane. Returns solution in terms of
    // line parameter t.
    std::vector<double> intersect(double3 origin, double3 direction) {
        double denominator = dot(normal, direction);
        if (denominator == 0) return {INFINITY};  // Plane is parallel to ray

        double t = -(distance + dot(normal, origin)) / denominator;
        if (t < 0) return {INFINITY}; // Triangle is 'behind' the ray
        return {t};
    }

    double3 get_normal_of(double3 point) {
        return normal;
    }
};

// Compute on which side a point lies, relative to the line defined by two
// points. This is the sign that you would get by computing
// normal.dot(point) - distance.
double sign(double3 point, double3 a, double3 b, double3 normal) {
    double3 candidate = subtract(point, a);
    double3 edge = subtract(b, a);
    return dot(normal, cross_product(candidate, edge));
}

struct Triangle : Object {
    double3 a;
    double3 b;
    double3 c;
    Plane plane;

    Triangle(double3 v_a, double3 v_b, double3 v_c, const double3& v_color, double v_specular, double v_reflective)
    : Object(v_color, v_specular, v_reflective, TRIANGLE) {
        a = v_a;
        b = v_b;
        c = v_c;

        // Computes plane containing triangle. Note this is redundant
        // with vertex information but allows us to cache computation.
        double3 ab = subtract(b, a);
        double3 ac = subtract(c, a);
        double3 normal = cross_product(ab, ac);
        double magnitude = length(normal);
        if (magnitude == 0) {
            std::cerr << "Error: Triangle is degenerate (i.e., points are collinear)" << std::endl;
            return;
        }
        normal = multiply(1.0f / magnitude, normal);

        double distance = -dot(normal, a);
        plane = Plane(normal, distance, v_color, v_specular, v_reflective);
    }

    // Computes intersection of a ray with triangle. Returns solution in terms
    // of line parameter t.
    std::vector<double> intersect(double3 origin, double3 direction) {
        // Find intersection between ray and plane
        double t = plane.intersect(origin, direction)[0];
        double3 point = add(origin, multiply(t, direction));

        // Check if the intersection lies in the triangle. NOTE the ordering of
        // points must be consistently clockwise OR counter clockwise, in this
        // calculation.
        double sign0 = sign(point, a, b, plane.normal);
        double sign1 = sign(point, b, c, plane.normal);
        double sign2 = sign(point, c, a, plane.normal);

        bool is_inside_cw = (sign0 > 0) && (sign1 > 0) && (sign2 > 0);
        bool is_inside_ccw = (sign0 < 0) && (sign1 < 0) && (sign2 < 0);
        if (is_inside_cw || is_inside_ccw) {
            return {t};
        }

        // If not in the triangle, ray does not intersect
        return {INFINITY};
    }

    double3 get_normal_of(double3 point) {
        return plane.get_normal_of(point);
    }
};

enum LightType {AMBIENT, POINT, DIRECTIONAL};

struct Light {
    LightType ltype;
    double intensity;
    double3 position;

    Light() {}

    Light(LightType v_ltype, double v_intensity, const double3& v_position) {
        ltype = v_ltype;
        intensity = v_intensity;
        position = v_position;
    }
};

struct Camera {
    double3 position;
    double33 rotation;

    Camera(double3 v_position, double33 v_rotation) {
        position = v_position;
        rotation = v_rotation;
    }
};

struct Scene {
    std::vector<Object*> objects;
    std::vector<Light> lights;

    Scene(std::vector<Object*> v_objects, std::vector<Light> v_lights) {
        objects = v_objects;
        lights = v_lights;
    }
};

// Convert 2d pixel coordinates to 3d viewport coordinates.
double3 canvas_to_viewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (double) x / width, (double) y / height, 1 };
}

// Holds intersection information for a raycast
struct Intersection {
    Object object;
    double3 point;
    double3 normal;
    bool is_valid;

    Intersection() {
        is_valid = false;
    }

    Intersection(Object v_object, double3 v_point, double3 v_normal) {
        object = v_object;
        point = v_point;
        normal = v_normal;
        is_valid = true;
    }
};

// Find the closest intersection between a ray and the spheres in the scene.
Intersection closest_intersection(
    double3 origin,
    double3 direction,
    double min_t,
    double max_t,
    Scene scene
) {
    double closest_t = INFINITY;
    Object* closest_object;

    for (int i = 0; i < scene.objects.size(); i++) {
        Object* object = scene.objects[i];

        std::vector<double> ts = object->intersect(origin, direction);
        for (int j = 0; j < ts.size(); j++) {
            if (ts[j] < closest_t && min_t < ts[j] && ts[j] < max_t) {
                closest_t = ts[j];
                closest_object = object;
            }
        }
    }

    // sets intersection.is_valid=false
    if (closest_t == INFINITY) return Intersection();

    // sets intersection.is_valid = true
    double3 point = add(origin, multiply(closest_t, direction));
    return Intersection(*closest_object, point, closest_object->get_normal_of(point));
}

// Compute the reflection of a ray on a surface defined by its normal
double3 reflect_ray(double3 ray, double3 normal) {
    return subtract(multiply(2 * dot(ray, normal), normal), ray);
}

// Compute lighting for the scene
double compute_lighting(double3 point, double3 normal, double3 view, double specular, Scene scene) {
    double intensity = 0;
    if (abs(length(normal) - 1.0f) > 0.0001f) {
        std::cerr << "Error: Normal is not length 1 (" << length(normal) << ")" << std::endl;
        return INFINITY;
    }

    double length_v = length(view);

    for (int i = 0; i < scene.lights.size(); i++) {
        Light light = scene.lights[i];
        if (light.ltype == AMBIENT) {
            intensity += light.intensity;
        } else {
            double3 vec_l;
            double shadow_t_max;

            if (light.ltype == POINT) {
                vec_l = subtract(light.position, point);
                shadow_t_max = 1.0f;
            } else {  // Light.DIRECTIONAL
                vec_l = light.position;
                shadow_t_max = INFINITY;
            }

            // Shadow check
            Intersection intersection = closest_intersection(point, vec_l, EPSILON, shadow_t_max, scene);
            if (intersection.is_valid) continue;

            // Diffuse
            double n_dot_l = dot(normal, vec_l);
            if (n_dot_l > 0) {
                intensity += light.intensity * n_dot_l / (length(vec_l));
            }

            // Specular, where vec_r is the 'perfect' reflection ray
            if (specular != -1) {
                double3 vec_r = reflect_ray(vec_l, normal);
                double r_dot_v = dot(vec_r, view);
                if (r_dot_v > 0) {
                    intensity += light.intensity * pow(r_dot_v / (length(vec_r) * length_v), specular);
                }
            }
        }
    }

    return intensity;
}

// Traces a ray against the spheres in the scene
double3 trace_ray(
    double3 origin,
    double3 direction,
    double min_t,
    double max_t,
    int32_t recursion_depth,
    Scene scene
) {
    Intersection intersection = closest_intersection(origin, direction, min_t, max_t, scene);
    if (!intersection.is_valid) {
        return BACKGROUND_COLOR;
    }

    double intensity = compute_lighting(
        intersection.point,
        intersection.normal,
        multiply(-1, direction),
        intersection.object.specular,
        scene);
    double3 local_color = multiply(intensity, intersection.object.color);

    // If we hit the recursion limit or the sphere is not reflective, finish
    double reflective = intersection.object.reflective;
    if (recursion_depth <= 0 || reflective <= 0.0f) {
        return local_color;
    }

    // Compute the reflected color
    double3 reflected_ray = reflect_ray(multiply(-1.0f, direction), intersection.normal);
    const double3 reflected_color = trace_ray(intersection.point, reflected_ray, EPSILON, INFINITY, recursion_depth - 1, scene);

    return add(multiply(1.0f - reflective, local_color), multiply(reflective, reflected_color));
}


int32_t main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    // Define camera settings
    double3 position = {3, 0, 1};
    double33 rotation = {{
        {{0.7071, 0, -0.7071}},
        {{     0, 1,       0}},
        {{0.7071, 0,  0.7071}}
    }};
    Camera camera = Camera(position, rotation);

    // Define scene
    std::vector<Object*> objects = {
        new Sphere({0, -1.0f, 3.0f}, 1.0f, {255, 0, 0}, 500.0f, 0.2f),
        new Sphere({-2.0f, 0, 4.0f}, 1.0f, {0, 255, 0}, 10.0f, 0.4f),
        new Sphere({2.0f, 0, 4.0f}, 1.0f, {0, 0, 255}, 500.0f, 0.3f),
        new Sphere({0, -5001.0f, 0}, 5000.0f, {255, 255, 0}, 1000.0f, 0.5f),
        new Triangle({1.0, 0.0, 5.0}, {-1.0, 0.0, 5.0}, {0.0, 2.0, 4.0}, {0, 255, 255}, 500.0, 0.4),
        new Triangle({1.0, 0.0, 5.0}, {0.0, 2.0, 4.0}, {0.0, 2.0, 6.0}, {0, 255, 255}, 500.0, 0.4),
        new Triangle({1.0, 0.0, 5.0}, {-1.0, 0.0, 5.0}, {0.0, 2.0, 6.0}, {0, 255, 255}, 500.0, 0.4),
        new Triangle({-1.0, 0.0, 5.0}, {0.0, 2.0, 4.0}, {0.0, 2.0, 6.0}, {0, 255, 255}, 500.0, 0.4)
    };

    std::vector<Light> lights = {
        Light(AMBIENT, 0.2f, {0, 0, 0}),
        Light(POINT, 0.6f, {2, 1, 0}),
        Light(DIRECTIONAL, 0.2f, {1, 4, 4})
    };
    Scene scene = Scene(objects, lights);

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            double3 direction = matmul(camera.rotation, canvas_to_viewport(x, y, width, height));
            double3 color = trace_ray(camera.position, direction, 1.0f, INFINITY, 3, scene);
            put_pixel(data, width, height, x, y, clamp(color));
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}