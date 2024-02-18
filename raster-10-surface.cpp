/*
Raster 10
=========
Fill triangles, and remove hidden surfaces

```bash
g++ raster-10-surface.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raster-10.html
*/
#include "bmp.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;
const int32_t WIDTH = 600;
const int32_t HEIGHT = 600;
const int VIEWPORT_SIZE = 1;
const int PROJECTION_PLANE_Z = 1;

// Canvas

bool PutPixel(
    uint8_t data[WIDTH * HEIGHT][3],
    int32_t x,
    int32_t y,
    const rgb color
) {
    // Translate [-W/2, W/2] into [0, W]. Ditto for H
    x = WIDTH / 2 + x;
    y = HEIGHT / 2 - y - 1;

    // Checks that the pixel is in bound
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) {
        // std::cerr << "Error: Attempted to write out-of-bounds pixel (" << x << ", " << y << ")." << std::endl;
        return false;
    }

    // Write into data, which is a flattened array where width * height in 1d
    int32_t offset = y * WIDTH + x;
    data[offset][0] = color[0];
    data[offset][1] = color[1];
    data[offset][2] = color[2];
    return true;
}

void Clear(uint8_t data[][3]) {
    for (int x = 0; x < WIDTH; x++) {
        for (int y = 0; y < HEIGHT; y++) {
            for (int i = 0; i < 3; i++) {
                data[y * WIDTH + x][i] = 255;
            }
        }
    }
}

// Depth buffer

bool UpdateDepthBufferIfCloser(int x, int y, float inv_z, float depth_buffer[WIDTH * HEIGHT]) {
    x = std::round(((float) WIDTH)/2 + x);
    y = std::round(((float) HEIGHT)/2 - y - 1);

    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) {
        return false;
    }

    int offset = x + WIDTH*y;
    if (depth_buffer[offset] == 0 || depth_buffer[offset] < inv_z) {
        depth_buffer[offset] = inv_z;
        return true;
    }
    return false;
}

// Data model

// A Point.
struct Point {
    int x;
    int y;

    Point(const int v_x, const int v_y) {
        x = v_x;
        y = v_y;
    }
};

// A 3D vertex.
struct Vertex {
    float x;
    float y;
    float z;

    Vertex() {}

    Vertex(const float v_x, const float v_y, const float v_z) {
        x = v_x;
        y = v_y;
        z = v_z;
    }
};

struct Vertex4 {
    float x;
    float y;
    float z;
    float w;

    Vertex4() {}

    Vertex4(Vertex vertex) {
        x = vertex.x;
        y = vertex.y;
        z = vertex.z;
        w = 1;
    }

    Vertex4(const float v_x, const float v_y, const float v_z, const float v_w = 1) {
        x = v_x;
        y = v_y;
        z = v_z;
        w = v_w;
    }
};

struct Mat4x4 {
    std::array<std::array<float, 4>, 4> data;

    Mat4x4() {}

    Mat4x4(const std::array<std::array<float, 4>, 4> v_data) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                data[i][j] = v_data[i][j];
            }
        }
    }
};

const Mat4x4 Identity4x4 = Mat4x4({{
    {{1, 0, 0, 0}},
    {{0, 1, 0, 0}},
    {{0, 0, 1, 0}},
    {{0, 0, 0, 1}}
}});

struct Triangle {
    std::array<int, 3> indexes;
    rgb color;

    Triangle(std::array<int, 3> v_indexes, rgb v_color) {
        indexes = v_indexes;
        color = v_color;
    }
};

struct Model {
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    Vertex bounds_center;
    float bounds_radius;
    bool is_valid;

    Model() {
        is_valid = false;
    }

    Model(const std::vector<Vertex>& v_vertices, const std::vector<Triangle>& v_triangles, Vertex v_bounds_center, float v_bounds_radius) {
        vertices = v_vertices;
        triangles = v_triangles;
        bounds_center = v_bounds_center;
        bounds_radius = v_bounds_radius;
        is_valid = true;
    }
};

struct Instance {
    Model model;
    Mat4x4 transform;
    float scale;

    Instance(Model v_model, Mat4x4 v_transform, float v_scale = 1) {
        model = v_model;
        transform = v_transform;
        scale = v_scale;
    }
};

// A clipping plane
struct Plane {
    Vertex normal;
    float distance;

    Plane(Vertex v_normal, float v_distance) {
        normal = v_normal;
        distance = v_distance;
    }
};

struct Camera {
    Vertex position;
    Mat4x4 orientation;
    std::vector<Plane> clipping_planes;

    Camera(Vertex v_position, Mat4x4 v_orientation) {
        position = v_position;
        orientation = v_orientation;
    }
};

// Linear algebra

Vertex Multiply(float k, Vertex vec) {
    return Vertex(k * vec.x, k * vec.y, k * vec.z);
}

// Note that the Javascript version of Dot also ignores the final w coordinate.
// This assumes that the coordinates have been converted into homoegenous
// coordinates, so the last v2.w is 1.
float Dot(Vertex v1, Vertex4 v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

float Dot(Vertex v1, Vertex v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

Vertex Add(Vertex v1, Vertex v2) {
    return Vertex(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vertex Cross(Vertex v1, Vertex v2) {
    return Vertex(
        v1.y*v2.z - v1.z*v2.y,
        v1.z*v2.x - v1.x*v2.z,
        v1.x*v2.y - v1.y*v2.x);
}

Mat4x4 MakeOYRotationMatrix(float degrees) {
    float cosine = cos(degrees * M_PI / 180);
    float sine = sin(degrees * M_PI / 180);

    return Mat4x4({{
        {{cosine, 0, -sine, 0}},
        {{0, 1, 0, 0}},
        {{sine, 0, cosine, 0}},
        {{0, 0, 0, 1}}
    }});
}

Mat4x4 MakeTranslationMatrix(Vertex translation) {
    return Mat4x4({{
        {{1, 0, 0, translation.x}},
        {{0, 1, 0, translation.y}},
        {{0, 0, 1, translation.z}},
        {{0, 0, 0, 1}}
    }});
}

Mat4x4 MakeScalingMatrix(float scale) {
    return Mat4x4({{
        {{scale, 0, 0, 0}},
        {{0, scale, 0, 0}},
        {{0, 0, scale, 0}},
        {{0, 0, 0, 1}}
    }});
}

Vertex4 MultiplyMV(Mat4x4 mat4x4, Vertex4 vec4) {
    std::array<float, 4> result = {0, 0, 0, 0};
    std::array<float, 4> vec = {vec4.x, vec4.y, vec4.z, vec4.w};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i] += mat4x4.data[i][j]*vec[j];
        }
    }

    return Vertex4(result[0], result[1], result[2], result[3]);
}

Vertex MultiplyMV(Mat4x4 mat4x4, Vertex vec3) {
    Vertex4 out = MultiplyMV(mat4x4, Vertex4(vec3));
    return Vertex(out.x, out.y, out.z);
}

Mat4x4 MultiplyMM4(Mat4x4 matA, Mat4x4 matB) {
  Mat4x4 result = Mat4x4({{
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}}
  }});

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        result.data[i][j] += matA.data[i][k]*matB.data[k][j];
      }
    }
  }

  return result;
}

Mat4x4 Transposed(Mat4x4 mat) {
  Mat4x4 result;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      result.data[i][j] = mat.data[j][i];
    }
  }
  return result;
}

template <typename T>
void Shuffle(std::vector<T> vec) {
    // initialize random number generator, uniform probability over all items
    std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<std::size_t> distribution(0, vec.size() - 1);

    for (int i = vec.size() - 1; i > 0; --i) {
        std::size_t number = distribution(generator);
        std::swap(vec[i], vec[rand]);
    }
}

rgb MultiplyColor(float k, rgb color) {
    return {
        (uint8_t) std::round(k * (float)color[0]),
        (uint8_t) std::round(k * (float)color[1]),
        (uint8_t) std::round(k * (float)color[2])
    };
}

Mat4x4 BuildTransformMatrix(Vertex v_position, Mat4x4 v_orientation = Identity4x4, float v_scale = 1) {
    return MultiplyMM4(MakeTranslationMatrix(v_position), MultiplyMM4(v_orientation, MakeScalingMatrix(v_scale)));
}

// Rasterization code

std::vector<float> Interpolate(int i0, float d0, int i1, float d1) {
    if (i0 == i1) {
        return {(float) d0};
    }

    std::vector<float> values = {};
    float a = ((float) d1 - d0) / (i1 - i0);
    float d = d0;
    for (int i = i0; i <= i1; i++) {
        values.push_back(d);
        d += a;
    }

    return values;
}

std::vector<int> FloatVectorToInt(const std::vector<float>& float_vec) {
    std::vector<int> int_vec(float_vec.size());
    for (size_t i = 0; i < float_vec.size(); ++i) {
        int_vec[i] = static_cast<int>(std::round(float_vec[i]));
    }
    return int_vec;
}

std::vector<int> Interpolate(int i0, int d0, int i1, int d1) {
    return FloatVectorToInt(Interpolate(i0, (float) d0, i1, (float) d1));
}

void DrawLine(uint8_t data[][3], Point p0, Point p1, rgb color) {
    int dx = p1.x - p0.x;
    int dy = p1.y - p0.y;

    if (std::abs(dx) > std::abs(dy)) {
        // The line is horizontal-ish. Make sure it's left to right.
        if (dx < 0) std::swap(p0, p1);

        // Compute the Y values and draw.
        std::vector<int> ys = Interpolate(p0.x, p0.y, p1.x, p1.y);
        for (int x = p0.x; x <= p1.x; x++) {
            PutPixel(data, x, ys[x - p0.x], color);
        }
    } else {
        // The line is verical-ish. Make sure it's bottom to top.
        if (dy < 0) std::swap(p0, p1);

        // Compute the X values and draw.
        std::vector<int> xs = Interpolate(p0.y, p0.x, p1.y, p1.x);
        for (int y = p0.y; y <= p1.y; y++) {
            PutPixel(data, xs[y - p0.y], y, color);
        }
    }
}

// Converts 2D viewport coordinates to 2D canvas coordinates.
Point ViewportToCanvas(float vx, float vy) {
    return Point(
        (int) std::round(vx * WIDTH / VIEWPORT_SIZE),
        (int) std::round(vy * HEIGHT / VIEWPORT_SIZE));
}

Point ProjectVertex(Vertex v) {
    return ViewportToCanvas(
        v.x * PROJECTION_PLANE_Z / v.z,
        v.y * PROJECTION_PLANE_Z / v.z);
}

rgb SortedVertexIndexes(std::array<int, 3> vertex_indexes, std::vector<Point> projected) {
    rgb indexes = {0, 1, 2};

    if (projected[vertex_indexes[indexes[1]]].y < projected[vertex_indexes[indexes[0]]].y) { std::swap(indexes[0], indexes[1]); }
    if (projected[vertex_indexes[indexes[2]]].y < projected[vertex_indexes[indexes[0]]].y) { std::swap(indexes[0], indexes[2]); }
    if (projected[vertex_indexes[indexes[2]]].y < projected[vertex_indexes[indexes[1]]].y) { std::swap(indexes[1], indexes[2]); }

    return indexes;
}

Vertex ComputeTriangleNormal(Vertex v0, Vertex v1, Vertex v2) {
    Vertex v0v1 = Add(v1, Multiply(-1, v0));
    Vertex v0v2 = Add(v2, Multiply(-1, v0));
    return Cross(v0v1, v0v2);
}

template <typename T>
std::vector<T> Concat(std::vector<T> v1, std::vector<T> v2, int start = 0) {
    std::vector<T> out;
    out.reserve(v1.size() + v2.size() - 1); // pre-allocate
    std::copy(v1.begin() + start, v1.end(), std::back_inserter(out)); // ignore first element
    std::copy(v2.begin(), v2.end(), std::back_inserter(out)); // concat second array
    return out;
}

template <typename T>
std::tuple<std::vector<T>, std::vector<T>> EdgeInterpolate(int y0, T v0, int y1, T v1, int y2, T v2) {
    std::vector<T> v01 = Interpolate(y0, v0, y1, v1);
    std::vector<T> v12 = Interpolate(y1, v1, y2, v2);
    std::vector<T> v02 = Interpolate(y0, v0, y2, v2);
    return std::make_tuple(v02, Concat(v01, v12));
}

void RenderTriangle(uint8_t data[WIDTH*HEIGHT][3], float depth_buffer[WIDTH*HEIGHT], Triangle triangle, std::vector<Vertex> vertices, std::vector<Point> projected) {
    // Sort by projected point Y.
    rgb indexes = SortedVertexIndexes(triangle.indexes, projected);
    int i0 = indexes[0];
    int i1 = indexes[1];
    int i2 = indexes[2];

    Vertex v0 = vertices[triangle.indexes[i0]];
    Vertex v1 = vertices[triangle.indexes[i1]];
    Vertex v2 = vertices[triangle.indexes[i2]];

    // Compute triangle normal. Use the unsorted vertices, otherwise the winding of the points may change.
    Vertex normal = ComputeTriangleNormal(vertices[triangle.indexes[0]], vertices[triangle.indexes[1]], vertices[triangle.indexes[2]]);

    // Backface culling.
    Vertex vertex_to_camera = Multiply(-1, vertices[triangle.indexes[0]]);  // Should be Subtract(camera.position, vertices[triangle.indexes[0]])
    if (Dot(vertex_to_camera, normal) <= 0) {
        return;
    }

    // Get attribute values (X, 1/Z) at the vertices.
    Point p0 = projected[triangle.indexes[i0]];
    Point p1 = projected[triangle.indexes[i1]];
    Point p2 = projected[triangle.indexes[i2]];

    // Compute attribute values at the edges.
    std::tuple xout = EdgeInterpolate(p0.y, p0.x, p1.y, p1.x, p2.y, p2.x);
    std::vector<int> x02 = std::get<0>(xout);
    std::vector<int> x012 = std::get<1>(xout);
    std::tuple izout = EdgeInterpolate(p0.y, (float) 1.0/v0.z, p1.y, (float) 1.0/v1.z, p2.y, (float) 1.0/v2.z);
    std::vector<float> iz02 = std::get<0>(izout);
    std::vector<float> iz012 = std::get<1>(izout);


    // Determine which is left and which is right.
    std::vector<int> x_left;
    std::vector<int> x_right;
    std::vector<float> iz_left;
    std::vector<float> iz_right;
    int m = (int) std::floor(x02.size() / 2);
    if (x02[m] < x012[m]) {
        x_left = x02;
        x_right = x012;
        iz_left = iz02;
        iz_right = iz012;
    } else {
        x_left = x012;
        x_right = x02;
        iz_left = iz012;
        iz_right = iz02;
    }

    // Draw horizontal segments.
    for (int y = p0.y; y <= p2.y; y++) {
        int xl = x_left[y - p0.y];
        int xr = x_right[y - p0.y];

        // Interpolate attributes for this scanline.
        float zl = iz_left[y - p0.y];
        float zr = iz_right[y - p0.y];
        std::vector<float> zscan = Interpolate(xl, zl, xr, zr);

        for (int x = xl; x <= xr; x++) {
            if (UpdateDepthBufferIfCloser(x, y, zscan[x - xl], depth_buffer)) {
                PutPixel(data, x, y, triangle.color);
            }
        }
    }

    rgb outline_color = MultiplyColor(0.75, triangle.color);
    DrawLine(data, p0, p1, outline_color);
    DrawLine(data, p1, p2, outline_color);
    DrawLine(data, p0, p2, outline_color);
}

// Clips a triangle against a plane. Adds output to triangles and vertices.
void ClipTriangle(Triangle triangle, Plane plane, std::vector<Triangle>& triangles, std::vector<Vertex> vertices) {
    Vertex v0 = vertices[triangle.indexes[0]];
    Vertex v1 = vertices[triangle.indexes[1]];
    Vertex v2 = vertices[triangle.indexes[2]];

    bool in0 = (Dot(plane.normal, v0) + plane.distance) > 0;
    bool in1 = (Dot(plane.normal, v1) + plane.distance) > 0;
    bool in2 = (Dot(plane.normal, v2) + plane.distance) > 0;

    int in_count = in0 + in1 + in2;
    if (in_count == 0) {
        // Nothing to do - the triangle is fully clipped out.
    } else if (in_count == 3) {
        // The triangle is fully in front of the plane.
        triangles.push_back(triangle);
    } else if (in_count == 2) {
        // The triangle has one vertex in. Output is one clipped triangle.
    } else if (in_count == 1) {
        // The triangle has two vertices in. Output is two clipped triangles.
    }
}

Model TransformAndClip(std::vector<Plane> clipping_planes, Model model, float scale, Mat4x4 transform) {
    // Transform the bounding sphere, and attempt early discard.
    Vertex center = MultiplyMV(transform, model.bounds_center);
    float radius = model.bounds_radius*scale;
    for (int p = 0; p < clipping_planes.size(); p++) {
        float distance = Dot(clipping_planes[p].normal, center) + clipping_planes[p].distance;
        if (distance < -radius) {
            return Model();
        }
    }

    // Apply modelview transform.
    std::vector<Vertex> vertices;
    for (int i = 0; i < model.vertices.size(); i++) {
        vertices.push_back(MultiplyMV(transform, model.vertices[i]));
        Vertex vertex = vertices[vertices.size() - 1];
    }

    // Clip the entire model against each successive plane.
    std::vector<Triangle> triangles(model.triangles); // performs copy
    for (int p = 0; p < clipping_planes.size(); p++) {
        std::vector<Triangle> new_triangles;
        for (int i = 0; i < triangles.size(); i++) {
            ClipTriangle(triangles[i], clipping_planes[p], new_triangles, vertices);
        }
        triangles = new_triangles;
    }

    return Model(vertices, triangles, center, model.bounds_radius);
}

void RenderModel(uint8_t data[WIDTH*HEIGHT][3], float depth_buffer[WIDTH*HEIGHT], Model model) {
    std::vector<Point> projected;
    for (int i = 0; i < model.vertices.size(); i++) {
        projected.push_back(ProjectVertex(model.vertices[i]));
    }
    for (int i = 0; i < model.triangles.size(); i++) {
        RenderTriangle(data, depth_buffer, model.triangles[i], model.vertices, projected);
    }
}

void RenderScene(uint8_t data[WIDTH*HEIGHT][3], float depth_buffer[WIDTH*HEIGHT], Camera camera, std::vector<Instance> instances) {
    Mat4x4 cameraMatrix = MultiplyMM4(Transposed(camera.orientation), MakeTranslationMatrix(Multiply(-1, camera.position)));

    for (int i = 0; i < instances.size(); i++) {
        Mat4x4 transform = MultiplyMM4(cameraMatrix, instances[i].transform);
        Model clipped = TransformAndClip(camera.clipping_planes, instances[i].model, instances[i].scale, transform);
        if (clipped.is_valid) RenderModel(data, depth_buffer, clipped);
    }
}

int main() {
    uint8_t data[WIDTH * HEIGHT][3];
    Clear(data);

    std::vector<Vertex> vertices = {
        Vertex(1, 1, 1),
        Vertex(-1, 1, 1),
        Vertex(-1, -1, 1),
        Vertex(1, -1, 1),
        Vertex(1, 1, -1),
        Vertex(-1, 1, -1),
        Vertex(-1, -1, -1),
        Vertex(1, -1, -1)
    };

    const rgb RED = {255, 0, 0};
    const rgb GREEN = {0, 255, 0};
    const rgb BLUE = {0, 0, 255};
    const rgb YELLOW = {255, 255, 0};
    const rgb PURPLE = {255, 0, 255};
    const rgb CYAN = {0, 255, 255};

    std::vector<Triangle> triangles = {
        Triangle({0, 1, 2}, RED),
        Triangle({0, 2, 3}, RED),
        Triangle({1, 5, 6}, YELLOW),
        Triangle({1, 6, 2}, YELLOW),
        Triangle({2, 6, 7}, CYAN),
        Triangle({2, 7, 3}, CYAN),
        Triangle({4, 0, 3}, GREEN),
        Triangle({4, 1, 0}, PURPLE),
        Triangle({4, 3, 7}, GREEN),
        Triangle({4, 5, 1}, PURPLE),
        Triangle({5, 4, 7}, BLUE),
        Triangle({5, 7, 6}, BLUE),
    };

    Model cube = Model(vertices, triangles, Vertex(0, 0, 0), sqrt(3));

    std::vector<Instance> instances = {
        // TODO: scale is passed twice
        Instance(cube, BuildTransformMatrix(Vertex(-1.5, 0, 7), Identity4x4, 0.75), 0.75),
        Instance(cube, BuildTransformMatrix(Vertex(1.25, 2.5, 7.5), MakeOYRotationMatrix(195))),
        Instance(cube, BuildTransformMatrix(Vertex(0, 0, -10), MakeOYRotationMatrix(195))),
    };

    Camera camera = Camera(Vertex(-3, 1, 2), MakeOYRotationMatrix(-30));
    float s2 = 1 / sqrt(2);
    camera.clipping_planes = {
        Plane(Vertex(0, 0, 1), -1), // Near
        Plane(Vertex(s2, 0, s2), 0), // Left
        Plane(Vertex(-s2, 0, s2), 0), // Right
        Plane(Vertex(0, -s2, s2), 0), // Top
        Plane(Vertex(0, s2, s2), 0), // Bottom
    };

    float depth_buffer[WIDTH*HEIGHT];
    for (int i = 0; i < WIDTH*HEIGHT; i++) {
        depth_buffer[i] = 0;
    }

    RenderScene(data, depth_buffer, camera, instances);

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, WIDTH, HEIGHT)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}