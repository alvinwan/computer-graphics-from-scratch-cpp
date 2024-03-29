#include <iostream>
#include <fstream>
#include <vector>

#pragma pack(push, 1)  // Ensure no padding is added by the compiler

// BMP file header structure
struct BMPFileHeader {
    uint16_t file_type;  // File type always BM which is 0x4D42
    uint32_t file_size;  // Size of the file (in bytes)
    uint16_t reserved1;  // Reserved, always 0
    uint16_t reserved2;  // Reserved, always 0
    uint32_t offset;     // Start position of pixel data
};

// BMP information header structure
struct BMPInfoHeader {
    uint32_t size;              // Size of the header (in bytes)
    int32_t width;              // Image width (in pixels)
    int32_t height;             // Image height (in pixels)
    uint16_t planes;            // Number of color planes (always 1)
    uint16_t bit_count;         // Number of bits per pixel
    uint32_t compression;       // Compression type (0 for uncompressed)
    uint32_t size_image;        // Image size in bytes (0 for uncompressed)
    int32_t x_pixels_per_meter; // Horizontal resolution in pixels per meter
    int32_t y_pixels_per_meter; // Vertical resolution in pixels per meter
    uint32_t colors_used;       // Number of colors used in the palette
    uint32_t colors_important;  // Number of important colors (0 - all important)
};

#pragma pack(pop)

bool write_bmp_file(const std::string& filename, const uint8_t data[][3], int32_t width, int32_t height) {
    // Checks that the dimensions are 'safe' for BMP
    if (width % 4 != 0) {
        std::cerr << "Error: Array width (" << width << ") must be divisible by 4." << std::endl;
        return false;
    }

    // Open the file in binary mode
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'." << std::endl;
        return false;
    }

    BMPInfoHeader info_header;
    BMPFileHeader file_header;

    // Fill the BMPInfoHeader structure
    info_header.size = sizeof(BMPInfoHeader);
    info_header.width = width;
    info_header.height = height;
    info_header.planes = 1;
    info_header.bit_count = 24;   // Assuming 3 bytes per pixel (24-bit color)
    info_header.compression = 0; // No compression
    info_header.size_image = width * height * 3;   // Assuming 3 bytes per pixel (24-bit color)
    info_header.x_pixels_per_meter = 2835; // Default resolution
    info_header.y_pixels_per_meter = 2835;
    info_header.colors_used = 0;
    info_header.colors_important = 0;

    // Fill the BMPFileHeader structure
    file_header.file_type = 0x4D42;
    file_header.offset = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader);
    file_header.file_size = file_header.offset + info_header.size_image; // File size
    file_header.reserved1 = 0;
    file_header.reserved2 = 0;

    // Write the file headers
    file.write((const char*)&file_header, sizeof(file_header));
    file.write((const char*)&info_header, sizeof(info_header));

    // Write the pixel data
    for (int32_t y = height - 1; y >= 0; y--) {
        for (int32_t x = 0; x < width; x++) {
            const uint8_t* rgb = data[y * width + x];
            file.put(rgb[2]); // blue
            file.put(rgb[1]); // green
            file.put(rgb[0]); // red
        }
    }

    file.close();

    if (file.good()) {
        std::cout << "BMP file '" << filename << "' written successfully." << std::endl;
        return true;
    } else {
        std::cerr << "Error: Error writing data to file." << std::endl;
        return false;
    }
}

bool read_bmp_file(const std::string filename, uint8_t data[][3], int32_t width, int32_t height) {
    // Open the file in binary mode
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'." << std::endl;
        return false;
    }

    BMPFileHeader file_header;
    BMPInfoHeader info_header;

    // Read the file headers
    file.read((char*)&file_header, sizeof(file_header));
    file.read((char*)&info_header, sizeof(info_header));

    // Validate the BMP headers
    if (file_header.file_type != 0x4D42 || info_header.bit_count != 24) {
        std::cerr << "Error: Invalid BMP file." << std::endl;
        return false;
    }

    width = info_header.width;
    height = info_header.height;

    // Calculate padding for each row (must be divisible by 4)
    int padding = (4 - ((info_header.width * 3) % 4)) % 4;

    // Read the pixel data
    for (int32_t y = height - 1; y >= 0; y--) {
        for (int32_t x = 0; x < width; x++) {
            const uint8_t* rgb = data[y * width + x];
            file.read((char*)&rgb[2], 1); // blue
            file.read((char*)&rgb[1], 1); // green
            file.read((char*)&rgb[0], 1); // red
        }
        // Skip padding bytes
        file.seekg(padding, std::ios_base::cur);
    }

    file.close();

    if (file.good()) {
        std::cout << "BMP file '" << filename << "' read successfully." << std::endl;
        return true;
    } else {
        std::cerr << "Error: Error reading data from file." << std::endl;
        return false;
    }
}