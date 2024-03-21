/*
Raycast 08 - Constructive Solid Geometry
========================================
Adds support for constructive solid geometry. Note that porting the compile-time
polymorphism changes from Raytracer 07 actually slowed down this script. This
may be because of non-locality in the new Material struct.

Timing: 1.08s

```bash
g++ raytracer-08-csg.cpp -o main.out -std=c++20 -Ofast
OUT=true ./main.out
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

// Multiply 3x3 matrix and 3x1 vector
float3 MultiplyMV(float33 matrix, float3 vector) {
    float3 out = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i] += vector[j] * matrix[i][j];
        }
    }
    return out;
}

// Implements cross product
float3 CrossProduct(float3 v1, float3 v2) {
    return {
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0]
    };
}

rgb Clamp(float3 vec) {
    return {
        (uint8_t) std::round(std::clamp<float>(vec[0], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(vec[1], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(vec[2], 0, 255)),
    };
}

// Ray tracing

struct Object {
    float3 color;
    float specular;
    float reflective;

    Object() {}

    Object(const float3& v_color, float v_specular, float v_reflective) {
        color = v_color;
        specular = v_specular;
        reflective = v_reflective;
    }

    // Computes intersection of a ray with object. Returns solution in terms of
    // line parameter t.
    virtual std::vector<float> Intersect(float3 origin, float3 direction) {
        return {INFINITY};
    }

    // Returns normal for a given point on the surface of this object.
    virtual float3 GetNormalOf(float3 point) {
        return {0, 0, 0};
    }

    virtual bool HasPoint(float3 point) {
        return false;
    }
};

struct Sphere : Object {
    float3 center;
    float radius;

    Sphere() {}

    Sphere(const float3& v_center, float v_radius, const float3& v_color, float v_specular, float v_reflective)
    : Object(v_color, v_specular, v_reflective) {
        center = v_center;
        radius = v_radius;
    }

    // Computes intersection of a ray with sphere. Returns solution in terms of
    // line parameter t.
    std::vector<float> Intersect(float3 origin, float3 direction) {
        float3 difference = Subtract(origin, center);

        float a = DotProduct(direction, direction);
        float b = 2 * DotProduct(difference, direction);
        float c = DotProduct(difference, difference) - radius * radius;

        float discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            return {INFINITY, INFINITY};
        }

        float t0 = (-b + sqrt(discriminant)) / (2 * a);
        float t1 = (-b - sqrt(discriminant)) / (2 * a);
        return {std::min(t0, t1), std::max(t0, t1)};
    }

    float3 GetNormalOf(float3 point) {
        float3 normal = Subtract(point, center);
        normal = Multiply(1 / Length(normal), normal);
        return normal;
    }

    bool HasPoint(float3 point) {
        return std::abs(Length(Subtract(point, center)) - radius) < EPSILON;
    }
};

struct Plane : Object {
    float3 normal;
    float distance;

    Plane() {}

    Plane(float3 v_normal, float v_distance, const float3& v_color, float v_specular, float v_reflective)
    : Object(v_color, v_specular, v_reflective) {
        normal = v_normal;
        distance = v_distance;
    }

    // Computes intersection of a ray with plane. Returns solution in terms of
    // line parameter t.
    std::vector<float> Intersect(float3 origin, float3 direction) {
        float denominator = DotProduct(normal, direction);
        // Plane is parallel to ray (=0) or we're looking at back side (>0)
        if (denominator >= 0) return {INFINITY};

        float t = -(distance + DotProduct(normal, origin)) / denominator;
        if (t < 0) return {INFINITY}; // Triangle is 'behind' the ray
        return {t};
    }

    float3 GetNormalOf(float3 point) {
        return normal;
    }

    bool HasPoint(float3 point) {
        return std::abs(DotProduct(point, normal) - distance) < EPSILON;
    }
};

// Compute on which side a point lies, relative to the line defined by two
// points. This is the sign that you would get by computing
// normal.DotProduct(point) - distance.
float Sign(float3 point, float3 a, float3 b, float3 normal) {
    float3 candidate = Subtract(point, a);
    float3 edge = Subtract(b, a);
    return DotProduct(normal, CrossProduct(candidate, edge));
}

struct Triangle : Object {
    float3 a;
    float3 b;
    float3 c;
    Plane plane;

    Triangle(float3 v_a, float3 v_b, float3 v_c, const float3& v_color, float v_specular, float v_reflective)
    : Object(v_color, v_specular, v_reflective) {
        a = v_a;
        b = v_b;
        c = v_c;

        // Computes plane containing triangle. Note this is redundant
        // with vertex information but allows us to cache computation.
        float3 ab = Subtract(b, a);
        float3 ac = Subtract(c, a);
        float3 normal = CrossProduct(ab, ac);
        float magnitude = Length(normal);
        if (magnitude == 0) {
            std::cerr << "Error: Triangle is degenerate (i.e., points are collinear)" << std::endl;
            return;
        }
        normal = Multiply(1 / magnitude, normal);

        float distance = -DotProduct(normal, a);
        plane = Plane(normal, distance, v_color, v_specular, v_reflective);
    }

    // Computes intersection of a ray with triangle. Returns solution in terms
    // of line parameter t.
    std::vector<float> Intersect(float3 origin, float3 direction) {
        // Find intersection between ray and plane
        float t = plane.Intersect(origin, direction)[0];
        float3 point = Add(origin, Multiply(t, direction));

        // Check if the intersection lies in the triangle. NOTE the ordering of
        // points must be consistently clockwise OR counter clockwise, in this
        // calculation. If *any sign is negative, ray does not intersect the
        // triangle.
        if (Sign(point, a, b, plane.normal) >= 0) return {INFINITY};
        if (Sign(point, b, c, plane.normal) >= 0) return {INFINITY};
        if (Sign(point, c, a, plane.normal) >= 0) return {INFINITY};

        // The intersection point lies in the 'correct' half plane for every
        // edge, so it lies in the triangle.
        return {t};
    }


    float3 GetNormalOf(float3 point) {
        return plane.GetNormalOf(point);
    }

    bool HasPoint(float3 point) {
        // TODO: This implementation is incomplete. Repeats intersect above.
        return plane.HasPoint(point);
    }
};

enum Operation { AND, OR, MINUS };

struct CSG : Object {
    Object* object1;
    Object* object2;
    Operation operation;

    CSG() {}

    CSG(Object* v_object1, Object* v_object2, Operation v_operation, const float3& v_color, float v_specular, float v_reflective)
    : Object(v_color, v_specular, v_reflective) {
        object1 = v_object1;
        object2 = v_object2;
        operation = v_operation;
    }

    std::vector<float> Intersect(float3 origin, float3 direction) {
        std::vector<float> ts1 = object1->Intersect(origin, direction);
        std::vector<float> ts2 = object2->Intersect(origin, direction);

        // Objects without volume will only return a single intersection, such
        // as a plane or triangle. Objects *with volume will return an even
        // number of interesections.
        if (ts1.size() % 2 != 0 || ts2.size() % 2 != 0) {
            std::cerr << "Error: Both objects in a CSG must have volume" << std::endl;
            return {INFINITY, INFINITY};
        }

        // Each operation has a different response, depending on which of the
        // two objects the ray intersects.
        if (operation == AND) {
            // In an AND operation, *both volumes must be intersected
            if (ts1[0] != INFINITY && ts2[0] != INFINITY) {
                return {std::min(ts1[1], ts2[1]), std::max(ts1[0], ts2[0])};
            }
            return {INFINITY, INFINITY};
        } else if (operation == OR) {
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

    float3 GetNormalOf(float3 point) {
        if (object1->HasPoint(point)) {
            // No matter the operation, if the point of intersection lies on the
            // first object, always use the first object's normal, normally.
            return object1->GetNormalOf(point);
        }
        // Otherwise, the point lies on the second object. In this case, we need
        // to handle operation-by-operation.
        float3 normal = object2->GetNormalOf(point);
        if (operation == AND || operation == OR) {
            // If operation is AND or OR, use the normal, normally.
            return normal;
        }
        // If the operation is a MINUS, and the point of intersection lies on
        // the second object, the normal needs to be flipped so it's on the
        // inside of object2.
        return Multiply(-1, normal);
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

struct Camera {
    float3 position;
    float33 rotation;

    Camera(float3 v_position, float33 v_rotation) {
        position = v_position;
        rotation = v_rotation;
    }
};

struct Scene {
    std::vector<Object*> objects;
    std::vector<Light> lights;
    float3 background_color;

    Scene(std::vector<Object*> v_objects, std::vector<Light> v_lights, float3 v_background_color) {
        objects = v_objects;
        lights = v_lights;
        background_color = v_background_color;
    }
};

// Convert 2d pixel coordinates to 3d viewport coordinates.
float3 CanvasToViewport(int32_t x, int32_t y, int32_t width, int32_t height) {
    return { (float) x / width, (float) y / height, 1 };
}

// Holds intersection information for a raycast
struct Intersection {
    Object object;
    float3 point;
    float3 normal;
    bool is_valid;

    Intersection() {
        is_valid = false;
    }

    Intersection(Object v_object, float3 v_point, float3 v_normal) {
        object = v_object;
        point = v_point;
        normal = v_normal;
        is_valid = true;
    }
};

// Find the closest intersection between a ray and the spheres in the scene.
Intersection ClosestIntersection(
    float3 origin,
    float3 direction,
    float min_t,
    float max_t,
    Scene scene
) {
    float closest_t = INFINITY;
    Object* closest_object;

    for (int i = 0; i < scene.objects.size(); i++) {
        Object* object = scene.objects[i];

        std::vector<float> ts = object->Intersect(origin, direction);
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
    float3 point = Add(origin, Multiply(closest_t, direction));
    return Intersection(*closest_object, point, closest_object->GetNormalOf(point));
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
            Intersection intersection = ClosestIntersection(point, vec_l, EPSILON, shadow_t_max, scene);
            if (intersection.is_valid) continue;

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
    Intersection intersection = ClosestIntersection(origin, direction, min_t, max_t, scene);
    if (!intersection.is_valid) {
        return scene.background_color;
    }

    float intensity = ComputeLighting(
        intersection.point,
        intersection.normal,
        Multiply(-1, direction),
        intersection.object.specular,
        scene);
    float3 local_color = Multiply(intensity, intersection.object.color);

    // If we hit the recursion limit or the sphere is not reflective, finish
    float reflective = intersection.object.reflective;
    if (recursion_depth <= 0 || reflective <= 0) {
        return local_color;
    }

    // Compute the reflected color
    float3 reflected_ray = ReflectRay(Multiply(-1, direction), intersection.normal);
    // NOTE: Below addresses 'shadow' acne (i.e., flickering shadow), which you
    // can see here: https://imgur.com/a/ycB69zX. To fix this, move the starting
    // point of the reflection ray along the normal, to prevent self-
    // intersection. JS demos don't have this problem because JS uses FP64 by
    // default. Additionally, other spheres don't have this problem because only
    // the biggest radius=5000 one has reflection rays that are near parallel.
    float3 point = Add(intersection.point, Multiply(EPSILON, intersection.normal));
    const float3 reflected_color = TraceRay(point, reflected_ray, EPSILON, INFINITY, recursion_depth - 1, scene);

    return Add(Multiply(1 - reflective, local_color), Multiply(reflective, reflected_color));
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
    std::vector<Object*> objects = {
        new Sphere({0, -5001, 0}, 5000, {255, 255, 0}, 1000, 0.5),
        new CSG(
            new Sphere({-2, 0, 4}, 1, {0, 0, 0}, 0, 0),
            new Sphere({-2, 1, 4}, 1, {0, 0, 0}, 0, 0),
            OR, {0, 255, 0}, 10, 0.4
        ),
        new CSG(
            new Sphere({0, -1, 3}, 1, {0, 0, 0}, 0, 0),
            new Sphere({0, 0, 3}, 1, {0, 0, 0}, 0, 0),
            AND, {255, 0, 0}, 500, 0.2f
        ),
        new CSG(
            new Sphere({2, 0, 4}, 1, {0, 0, 0}, 0, 0),
            new Sphere({2, 1, 3}, 1, {0, 0, 0}, 0, 0),
            MINUS, {0, 0, 255}, 500, 0.3
        ),
    };

    std::vector<Light> lights = {
        Light(AMBIENT, 0.2, {0, 0, 0}),
        Light(POINT, 0.6, {2, 1, 0}),
        Light(DIRECTIONAL, 0.2, {1, 4, 4})
    };
    Scene scene = Scene(objects, lights, {0, 0, 0});

    for (int32_t x = -width / 2; x < width / 2; x++) {
        for (int32_t y = -height / 2; y < height / 2; y++)
        {
            float3 direction = MultiplyMV(camera.rotation, CanvasToViewport(x, y, width, height));
            float3 color = TraceRay(camera.position, direction, 1, INFINITY, 3, scene);
            PutPixel(data, width, height, x, y, Clamp(color));
        }
    }

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}