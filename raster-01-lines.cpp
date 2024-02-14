/*
Raster 01
=========
Implements a simpler rasterization program for lines. This is the initial, buggy
line drawing algorithm.

```bash
g++ raster-01-lines.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-01.html
*/
#include "bmp.h"

typedef std::array<float, 3> float3;
typedef std::array<uint8_t, 3> rgb;

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

// Rasterization code

void draw_line_broken(uint8_t data[][3], int width, int height, int x0, int y0, int x1, int y1, rgb color) {
    float a = (y1 - y0) / (x1 - x0);
    float y = y0;
    for (int x = x0; x <= x1; x++) {
        put_pixel(data, width, height, x, (int) std::round(y), color);
        y += a;
    }
}

void clear(uint8_t data[][3], int width, int height) {
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            for (int i = 0; i < 3; i++) {
                data[y * width + x][i] = 255;
            }
        }
    }
}

int main() {
    int32_t width = 600;
    int32_t height = 600;
    uint8_t data[width * height][3];

    clear(data, width, height);

    draw_line_broken(data, width, height, -200, -100, 240, 120, {0, 0, 0});
    draw_line_broken(data, width, height, -50, -200, 60, 240, {0, 0, 0});

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}