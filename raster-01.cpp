#include "bmp.h"

int32_t main() {
    // Sample usage: Write a 200x100 image from an array of ints
    std::vector<int32_t> data(200 * 100, 0x00FF00); // Green color (0, 255, 0)

    if (write_bmp_file("output.bmp", data, 200, 100)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}