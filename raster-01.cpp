// Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/02-basic-raytracing.html
#include "bmp.h"


bool put_pixel(
    uint8_t data[][3],
    int32_t width,
    int32_t height,
    int32_t x,
    int32_t y,
    uint8_t color[3]
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


int32_t main() {
    int32_t width = 8;
    int32_t height = 8;
    uint8_t data[width * height][3];
 
    // Fill up the viewport with white pixels. Test that our conversions are
    // correct, from viewport to canvas.
    uint8_t color[3] = {255, 255, 255};
    for (int x = -width / 2; x < width / 2; x++) {
        for (int y = -height / 2; y < height / 2; y++)
        {
            put_pixel(data, width, height, x, y, color);
        }
    }

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}