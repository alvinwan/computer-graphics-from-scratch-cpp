/*
Raster 04
=========
Adding triangle shading

```bash
g++ raster-04-shading.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raster-04.html
*/
#include "bmp.h"
#include <vector>
#include <algorithm>
#include <cmath>

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;
const int32_t WIDTH = 600;
const int32_t HEIGHT = 600;

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

struct Point {
    int x;
    int y;
    float h;

    Point(const int v_x, const int v_y, const float v_h) {
        x = v_x;
        y = v_y;
        h = v_h;
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

void draw_wireframe_triangle(uint8_t data[][3], Point p0, Point p1, Point p2, rgb color) {
    draw_line(data, p0, p1, color);
    draw_line(data, p1, p2, color);
    draw_line(data, p0, p2, color);
}

template <typename T>
std::vector<T> concat(std::vector<T> v1, std::vector<T> v2, int start = 0) {
    std::vector<T> out;
    out.reserve(v1.size() + v2.size() - 1); // pre-allocate
    std::copy(v1.begin() + start, v1.end(), std::back_inserter(out)); // ignore first element
    std::copy(v2.begin(), v2.end(), std::back_inserter(out)); // concat second array
    return out;
}

void draw_shaded_triangle(uint8_t data[][3], Point p0, Point p1, Point p2, rgb color) {
    if (p1.y < p0.y) std::swap(p0, p1);
    if (p2.y < p0.y) std::swap(p0, p2);
    if (p2.y < p1.y) std::swap(p1, p2);

    // Compute X coordinates and H values of the edges.
    std::vector<int> x01 = interpolate(p0.y, p0.x, p1.y, p1.x);
    std::vector<float> h01 = interpolate(p0.y, p0.h, p1.y, p1.h);

    std::vector<int> x12 = interpolate(p1.y, p1.x, p2.y, p2.x);
    std::vector<float> h12 = interpolate(p1.y, p1.h, p2.y, p2.h);

    std::vector<int> x02 = interpolate(p0.y, p0.x, p2.y, p2.x);
    std::vector<float> h02 = interpolate(p0.y, p0.h, p2.y, p2.h);

    // Merge the two short sides.
    std::vector<int> x012 = concat(x01, x12, 1);
    std::vector<float> h012 = concat(h01, h12, 1);

    // Determine which is left and which is right.
    std::vector<int> x_left;
    std::vector<int> x_right;
    std::vector<float> h_left;
    std::vector<float> h_right;
    int m = (int) std::floor(x02.size() / 2);
    if (x02[m] < x012[m]) {
        x_left = x02;
        x_right = x012;
        h_left = h02;
        h_right = h012;
    } else {
        x_left = x012;
        x_right = x02;
        h_left = h012;
        h_right = h02;
    }

    // Draw horizontal segments.
    for (int y = p0.y; y <= p2.y; y++) {
        int xl = x_left[y - p0.y] | 0;
        int xr = x_right[y - p0.y] | 0;
        std::vector<float> h_segment = interpolate(xl, h_left[y - p0.y], xr, h_right[y - p0.y]);

        for (int x = xl; x <= xr; x++) {

            put_pixel(data, x, y, multiply(h_segment[x - xl], color));
        }
    }
}

int main() {
    uint8_t data[WIDTH * HEIGHT][3];

    clear(data);

    Point p0 = Point(-200, -250, 0.3);
    Point p1 = Point(200, 50, 0.1);
    Point p2 = Point(20, 250, 1.0);

    draw_shaded_triangle(data, p0, p1, p2, {0, 255, 0});

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, WIDTH, HEIGHT)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}