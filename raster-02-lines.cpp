/*
Raster 02
=========
Fix rasterization of lines

```bash
g++ raster-02-lines.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raster-02.html
*/
#include "bmp.h"

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

    Point(const int v_x, const int v_y) {
        x = v_x;
        y = v_y;
    }
};

// Rasterization code

std::vector<int> interpolate(int i0, int d0, int i1, int d1) {
    if (i0 == i1) {
        return {d0};
    }

    std::vector<int> values = {};
    float a = ((float) d1 - d0) / (i1 - i0);
    float d = d0;
    for (int i = i0; i <= i1; i++) {
        values.push_back((int) std::round(d));
        d += a;
    }

    return values;
}

void draw_line(uint8_t data[WIDTH * HEIGHT][3], Point p0, Point p1, rgb color) {
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

int main() {
    uint8_t data[WIDTH * HEIGHT][3];

    clear(data);

    draw_line(data, Point(-200, -100), Point(240, 120), {0, 0, 0});
    draw_line(data, Point(-50, -200), Point(60, 240), {0, 0, 0});

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, WIDTH, HEIGHT)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}