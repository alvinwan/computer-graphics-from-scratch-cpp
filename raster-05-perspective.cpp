/*
Raster 05
=========
Perspective projection with a cube

```bash
g++ raster-05-perspective.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raster-05.html
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

bool put_pixel(
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

void clear(uint8_t data[][3]) {
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

    Vertex(const float v_x, const float v_y, const float v_z) {
        x = v_x;
        y = v_y;
        z = v_z;
    }
};

// Linear algebra

rgb multiply(float k, rgb vec) {
    return {
        (uint8_t) std::round(std::clamp<float>(k * vec[0], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(k * vec[1], 0, 255)),
        (uint8_t) std::round(std::clamp<float>(k * vec[2], 0, 255))
    };
}

// Rasterization code

std::vector<float> interpolate(int i0, float d0, int i1, float d1) {
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

std::vector<int> float_vector_to_int(const std::vector<float>& float_vec) {
    std::vector<int> int_vec(float_vec.size());
    for (size_t i = 0; i < float_vec.size(); ++i) {
        int_vec[i] = static_cast<int>(std::round(float_vec[i]));
    }
    return int_vec;
}

std::vector<int> interpolate(int i0, int d0, int i1, int d1) {
    return float_vector_to_int(interpolate(i0, (float) d0, i1, (float) d1));
}

void draw_line(uint8_t data[][3], Point p0, Point p1, rgb color) {
    int dx = p1.x - p0.x;
    int dy = p1.y - p0.y;

    if (std::abs(dx) > std::abs(dy)) {
        // The line is horizontal-ish. Make sure it's left to right.
        if (dx < 0) std::swap(p0, p1);

        // Compute the Y values and draw.
        std::vector<int> ys = interpolate(p0.x, p0.y, p1.x, p1.y);
        for (int x = p0.x; x <= p1.x; x++) {
            put_pixel(data, x, ys[x - p0.x], color);
        }
    } else {
        // The line is verical-ish. Make sure it's bottom to top.
        if (dy < 0) std::swap(p0, p1);

        // Compute the X values and draw.
        std::vector<int> xs = interpolate(p0.y, p0.x, p1.y, p1.x);
        for (int y = p0.y; y <= p1.y; y++) {
            put_pixel(data, xs[y - p0.y], y, color);
        }
    }
}

// Converts 2D viewport coordinates to 2D canvas coordinates.
Point viewport_to_canvas(float vx, float vy) {
    return Point(
        (int) std::round(vx * WIDTH / VIEWPORT_SIZE),
        (int) std::round(vy * HEIGHT / VIEWPORT_SIZE));
}


Point project_vertex(Vertex v) {
    return viewport_to_canvas(
        v.x * PROJECTION_PLANE_Z / v.z,
        v.y * PROJECTION_PLANE_Z / v.z);
}

int main() {
    uint8_t data[WIDTH * HEIGHT][3];
    clear(data);

    Vertex vA = Vertex(-2, -0.5, 5);
    Vertex vB = Vertex(-2, 0.5, 5);
    Vertex vC = Vertex(-1, 0.5, 5);
    Vertex vD = Vertex(-1, -0.5, 5);

    Vertex vAb = Vertex(-2, -0.5, 6);
    Vertex vBb = Vertex(-2, 0.5, 6);
    Vertex vCb = Vertex(-1, 0.5, 6);
    Vertex vDb = Vertex(-1, -0.5, 6);

    rgb RED = {255, 0, 0};
    rgb GREEN = {0, 255, 0};
    rgb BLUE = {0, 0, 255};

    draw_line(data, project_vertex(vA), project_vertex(vB), BLUE);
    draw_line(data, project_vertex(vB), project_vertex(vC), BLUE);
    draw_line(data, project_vertex(vC), project_vertex(vD), BLUE);
    draw_line(data, project_vertex(vD), project_vertex(vA), BLUE);

    draw_line(data, project_vertex(vAb), project_vertex(vBb), RED);
    draw_line(data, project_vertex(vBb), project_vertex(vCb), RED);
    draw_line(data, project_vertex(vCb), project_vertex(vDb), RED);
    draw_line(data, project_vertex(vDb), project_vertex(vAb), RED);

    draw_line(data, project_vertex(vA), project_vertex(vAb), GREEN);
    draw_line(data, project_vertex(vB), project_vertex(vBb), GREEN);
    draw_line(data, project_vertex(vC), project_vertex(vCb), GREEN);
    draw_line(data, project_vertex(vD), project_vertex(vDb), GREEN);

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, WIDTH, HEIGHT)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}