// Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/02-basic-raytracing.html
#include "bmp.h"


void put_pixel(
    std::vector<int32_t>& data,
    int32_t width,
    int32_t height,
    int32_t x,
    int32_t y,
    uint8_t color[3]
) {
    if (x < 0 || x >= width || y < 0 || y >= height) {
        return;
    }

    // Convert each color value to a byte, ensuring values are within 0-255 range
    uint8_t red = static_cast<uint8_t>(color[0]);
    uint8_t green = static_cast<uint8_t>(color[1]);
    uint8_t blue = static_cast<uint8_t>(color[2]);

    // Combine the shifted values using bitwise OR to create the final chat
    data[y * width + x] = (red << 16) | (green << 8) | blue;
}


int32_t main() {
    uint32_t width = 8;
    uint32_t height = 8;
    std::vector<int32_t> data(width * height, 0);

    uint8_t color[3] = {0, 255, 0};
    put_pixel(data, width, height, 0, 0, color);

    if (write_bmp_file("output.bmp", data, width, height)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}