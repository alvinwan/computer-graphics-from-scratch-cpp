/*
Raster 08
=========
Use homogenous coordinates and generalized camera.

```bash
g++ raster-08-homogenous.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raster-08.html
*/
#include "bmp.h"
#include <vector>
#include <algorithm>
#include <cmath>

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
        std::cerr << "Error: Attempted to write out-of-bounds pixel (" << x << ", " << y << ")." << std::endl;
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

    Vertex4(const float v_x, const float v_y, const float v_z, const float v_w) {
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
    int v0;
    int v1;
    int v2;
    rgb color;

    Triangle(int v_v0, int v_v1, int v_v2, rgb v_color) {
        v0 = v_v0;
        v1 = v_v1;
        v2 = v_v2;
        color = v_color;
    }
};

struct Model {
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;

    Model() {}

    Model(const std::vector<Vertex>& v_vertices, const std::vector<Triangle>& v_triangles) {
        vertices = v_vertices;
        triangles = v_triangles;
    }
};

struct Camera {
    Vertex position;
    Mat4x4 orientation;

    Camera(Vertex v_position, Mat4x4 v_orientation) {
        position = v_position;
        orientation = v_orientation;
    }
};

// Linear algebra

Vertex Multiply(float k, Vertex vec) {
    return Vertex(k * vec.x, k * vec.y, k * vec.z);
}

Vertex Add(Vertex v1, Vertex v2) {
    return Vertex(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
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
  Mat4x4 result = Mat4x4({{
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}},
    {{0, 0, 0, 0}}
  }});

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      result.data[i][j] = mat.data[j][i];
    }
  }
  return result;
}

// Rasterization code

// Technically, should create a header file so this struct can join the other
// structs above.
struct Instance {
    Model model;
    Vertex position;
    Mat4x4 orientation;
    float scale;
    Mat4x4 transform;

    Instance(Model v_model, Vertex v_position, Mat4x4 v_orientation = Identity4x4, float v_scale = 1) {
        model = v_model;
        position = v_position;
        orientation = v_orientation;
        scale = v_scale;
        transform = MultiplyMM4(MakeTranslationMatrix(v_position), MultiplyMM4(v_orientation, MakeScalingMatrix(v_scale)));
    }
};

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

void DrawWireframeTriangle(uint8_t data[WIDTH * HEIGHT][3], Point p0, Point p1, Point p2, rgb color) {
    DrawLine(data, p0, p1, color);
    DrawLine(data, p1, p2, color);
    DrawLine(data, p0, p2, color);
}

// Converts 2D viewport coordinates to 2D canvas coordinates.
Point ViewportToCanvas(float vx, float vy) {
    return Point(
        (int) std::round(vx * WIDTH / VIEWPORT_SIZE),
        (int) std::round(vy * HEIGHT / VIEWPORT_SIZE));
}

Point ProjectVertex(Vertex4 v) {
    return ViewportToCanvas(
        v.x * PROJECTION_PLANE_Z / v.z,
        v.y * PROJECTION_PLANE_Z / v.z);
}

void RenderTriangle(uint8_t data[WIDTH*HEIGHT][3], Triangle triangle, std::vector<Point> projected) {
    DrawWireframeTriangle(
        data,
        projected[triangle.v0],
        projected[triangle.v1],
        projected[triangle.v2],
        triangle.color);
}

void RenderModel(uint8_t data[WIDTH*HEIGHT][3], Model model, Mat4x4 transform) {
    std::vector<Point> projected;
    for (int i = 0; i < model.vertices.size(); i++) {
        Vertex vertex = model.vertices[i];
        Vertex4 vertexH = Vertex4(vertex.x, vertex.y, vertex.z, 1);
        projected.push_back(ProjectVertex(MultiplyMV(transform, vertexH)));
    }
    for (int i = 0; i < model.triangles.size(); i++) {
        RenderTriangle(data, model.triangles[i], projected);
    }
}

void RenderScene(uint8_t data[WIDTH*HEIGHT][3], Camera camera, std::vector<Instance> instances) {
    Mat4x4 cameraMatrix = MultiplyMM4(Transposed(camera.orientation), MakeTranslationMatrix(Multiply(-1, camera.position)));

    for (int i = 0; i < instances.size(); i++) {
        Mat4x4 transform = MultiplyMM4(cameraMatrix, instances[i].transform);
        RenderModel(data, instances[i].model, transform);
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
        Triangle(0, 1, 2, RED),
        Triangle(0, 2, 3, RED),
        Triangle(4, 0, 3, GREEN),
        Triangle(4, 3, 7, GREEN),
        Triangle(5, 4, 7, BLUE),
        Triangle(5, 7, 6, BLUE),
        Triangle(1, 5, 6, YELLOW),
        Triangle(1, 6, 2, YELLOW),
        Triangle(4, 5, 1, PURPLE),
        Triangle(4, 1, 0, PURPLE),
        Triangle(2, 6, 7, CYAN),
        Triangle(2, 7, 3, CYAN)
    };

    Camera camera = Camera(Vertex(-3, 1, 2), MakeOYRotationMatrix(-30));

    Model cube = Model(vertices, triangles);

    std::vector<Instance> instances = {
        Instance(cube, Vertex(-1.5, 0, 7), Identity4x4, 0.75),
        Instance(cube, Vertex(1.25, 2, 7.5), MakeOYRotationMatrix(195), 1)
    };

    RenderScene(data, camera, instances);

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, WIDTH, HEIGHT)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}