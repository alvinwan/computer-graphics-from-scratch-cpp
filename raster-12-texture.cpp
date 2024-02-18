/*
Raster 12
=========
Add textures

```bash
g++ raster-12-texture.cpp -o main.out -std=c++20 -Ofast
./main.out
open output.bmp
```

Implementation for https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raster-12.html
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

const bool USE_VERTEX_NORMALS = true;
enum Shading {FLAT, GOURAUD, PHONG};
const Shading SHADING_MODEL = PHONG;

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

    Point() {}

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
    // Unlike the JS demos, I moved Texture from the Triangle to the Instance
    std::array<int, 3> indexes;
    rgb color;
    std::array<Vertex, 3> normals;
    std::array<Point, 3> uvs;

    Triangle(std::array<int, 3> v_indexes, rgb v_color, std::array<Vertex, 3> v_normals, std::array<Point, 3> v_uvs) {
        indexes = v_indexes;
        color = v_color;
        normals = v_normals;
        uvs = v_uvs;
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

struct Texture {
    uint8_t data[WIDTH * HEIGHT][3];

    Texture() {}

    Texture(std::string filename) {
        read_bmp_file(filename, data, WIDTH, HEIGHT);
    }

    rgb getTexel(float u, float v) {
        int iu = (int) std::floor(u * WIDTH);
        int iv = (int) std::floor(v * HEIGHT);

        int offset = iv * WIDTH + iu;
        return {data[offset][0], data[offset][1], data[offset][2]};
    }
};

struct Instance {
    // Unlike the JS demos, I moved Texture from the Triangle to the Instance
    Model model;
    Mat4x4 transform;
    Mat4x4 orientation;
    Texture texture;
    float scale;

    Instance(Model v_model, Mat4x4 v_transform, Mat4x4 v_orientation,  Texture v_texture, float v_scale = 1) {
        model = v_model;
        transform = v_transform;
        orientation = v_orientation;
        scale = v_scale;
        texture = v_texture;
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

enum LightType {AMBIENT, POINT, DIRECTIONAL};

struct Light {
    LightType type;
    float intensity;
    Vertex vector;

    Light(LightType v_type, float v_intensity, Vertex v_vector) {
        type = v_type;
        intensity = v_intensity;
        vector = v_vector;
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

Vertex Add(Vertex4 v1, Vertex v2) {
    return Vertex(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vertex Cross(Vertex v1, Vertex v2) {
    return Vertex(
        v1.y*v2.z - v1.z*v2.y,
        v1.z*v2.x - v1.x*v2.z,
        v1.x*v2.y - v1.y*v2.x);
}

float Magnitude(Vertex v1) {
    return sqrt(Dot(v1, v1));
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
    Vertex4 out = MultiplyMV(mat4x4, Vertex4(vec3.x, vec3.y, vec3.z, 1));
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
uint8_t Clamp(T value) {
    return (uint8_t) std::round(std::clamp<T>(value, 0, 255));
}

// Adds two colors.
rgb AddColor(rgb c1, rgb c2) {
    return {
        Clamp((int) c1[0] + (int) c2[0]),
        Clamp((int) c1[1] + (int) c2[1]),
        Clamp((int) c1[2] + (int) c2[2])
    };
}

rgb MultiplyColor(float k, rgb color) {
    return {
        Clamp(k * (float)color[0]),
        Clamp(k * (float)color[1]),
        Clamp(k * (float)color[2])
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

Point CanvasToViewport(float vx, float vy) {
    return Point(
        (int) std::round(vx * VIEWPORT_SIZE / WIDTH),
        (int) std::round(vy * VIEWPORT_SIZE / HEIGHT)
    );
}

Point ProjectVertex(Vertex v) {
    return ViewportToCanvas(
        v.x * PROJECTION_PLANE_Z / v.z,
        v.y * PROJECTION_PLANE_Z / v.z);
}

Vertex UnProjectVertex(int x, int y, float inv_z) {
    float oz = 1.0 / inv_z;
    float ux = x * oz / PROJECTION_PLANE_Z;
    float uy = y * oz / PROJECTION_PLANE_Z;
    Point p2d = CanvasToViewport(ux, uy);
    return Vertex(p2d.x, p2d.y, oz);
}

// Sort the points from bottom to top.
// Technically, sort the indexes to the vertex indexes in the triangle from bottom to top.
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

float ComputeIllumination(Vertex vertex, Vertex normal, Camera camera, std::vector<Light> lights) {
    float illumination = 0;
    for (int l = 0; l < lights.size(); l++) {
        Light light = lights[l];
        if (light.type == AMBIENT) {
            illumination += light.intensity;
            continue;
        }

        Vertex vl;
        if (light.type == DIRECTIONAL) {
            Mat4x4 cameraMatrix = Transposed(camera.orientation);
            Vertex rotated_light = MultiplyMV(cameraMatrix, light.vector);
            vl = rotated_light;
        } else if (light.type == POINT) {
            Mat4x4 cameraMatrix = MultiplyMM4(Transposed(camera.orientation), MakeTranslationMatrix(Multiply(-1, camera.position)));
            Vertex transformed_light = MultiplyMV(cameraMatrix, light.vector);
            vl = Add(transformed_light, Multiply(-1, vertex)); // light.vector - vertex
        }

        // Diffuse component.
        float cos_alpha = Dot(vl, normal) / (Magnitude(vl) * Magnitude(normal));
        if (cos_alpha > 0) {
            illumination += cos_alpha * light.intensity;
        }

        // Specular component.
        Vertex reflected = Add(Multiply(2*Dot(normal, vl), normal), Multiply(-1, vl));
        Vertex view = Add(camera.position, Multiply(-1, vertex));

        float cos_beta = Dot(reflected, view) / (Magnitude(reflected) * Magnitude(view));
        if (cos_beta > 0) {
            float specular = 50;
            illumination += pow(cos_beta, specular) * light.intensity;
        }
    }
    return illumination;
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

void RenderTriangle(uint8_t data[WIDTH*HEIGHT][3], float depth_buffer[WIDTH*HEIGHT], Triangle triangle, std::vector<Vertex> vertices, std::vector<Point> projected, Camera camera, std::vector<Light> lights, Mat4x4 orientation, Texture& texture) {
    // Compute triangle normal. Use the unsorted vertices, otherwise the winding of the points may change.
    Vertex normal = ComputeTriangleNormal(vertices[triangle.indexes[0]], vertices[triangle.indexes[1]], vertices[triangle.indexes[2]]);

    // Backface culling.
    Vertex vertex_to_camera = Multiply(-1, vertices[triangle.indexes[0]]);  // Should be Subtract(camera.position, vertices[triangle.indexes[0]])
    if (Dot(vertex_to_camera, normal) <= 0) {
        return;
    }

    // Sort by projected point Y.
    rgb indexes = SortedVertexIndexes(triangle.indexes, projected);
    int i0 = indexes[0];
    int i1 = indexes[1];
    int i2 = indexes[2];

    Vertex v0 = vertices[triangle.indexes[i0]];
    Vertex v1 = vertices[triangle.indexes[i1]];
    Vertex v2 = vertices[triangle.indexes[i2]];

    // Get attribute values (X, 1/Z, U/Z, U/Z) at the vertices.
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

    // Apply texture -- assumes all triangles have a texture included
    std::tuple uz = EdgeInterpolate(
        p0.y, (float) triangle.uvs[i0].x / v0.z,
        p1.y, (float) triangle.uvs[i1].x / v1.z,
        p2.y, (float) triangle.uvs[i2].x / v2.z
    );
    std::vector<float> uz_left = std::get<0>(uz);
    std::vector<float> uz_right = std::get<1>(uz);
    std::tuple vz = EdgeInterpolate(
        p0.y, (float) triangle.uvs[i0].y / v0.z,
        p1.y, (float) triangle.uvs[i1].y / v1.z,
        p2.y, (float) triangle.uvs[i2].y / v2.z
    );
    std::vector<float> vz_left = std::get<0>(vz);
    std::vector<float> vz_right = std::get<1>(vz);

    // Use model normals instead of vertex normals
    Vertex normal0;
    Vertex normal1;
    Vertex normal2;
    if (USE_VERTEX_NORMALS) {
        Mat4x4 transform = MultiplyMM4(Transposed(camera.orientation), orientation);
        normal0 = MultiplyMV(transform, triangle.normals[i0]);
        normal1 = MultiplyMV(transform, triangle.normals[i1]);
        normal2 = MultiplyMV(transform, triangle.normals[i2]);
    } else {
        normal0 = normal;
        normal1 = normal;
        normal2 = normal;
    }

    std::vector<float> i02;
    std::vector<float> i012;
    std::vector<float> nx02;
    std::vector<float> nx012;
    std::vector<float> ny02;
    std::vector<float> ny012;
    std::vector<float> nz02;
    std::vector<float> nz012;
    float intensity;
    if (SHADING_MODEL == FLAT) {
        // Flat shading: compute lighting for the entire triangle
        Vertex center = Vertex((v0.x + v1.x + v2.x)/3.0, (v0.y + v1.y + v2.y)/3.0, (v0.z + v1.z + v2.z)/3.0);
        intensity = ComputeIllumination(center, normal0, camera, lights);
    } else if (SHADING_MODEL == GOURAUD) {
        // Gouraud shading: compute lighting at the vertices, and interpolate
        float i0 = ComputeIllumination(v0, normal0, camera, lights);
        float i1 = ComputeIllumination(v1, normal1, camera, lights);
        float i2 = ComputeIllumination(v2, normal2, camera, lights);
        std::tuple i0x = EdgeInterpolate(p0.y, i0, p1.y, i1, p2.y, i2);
        std::vector<float> i02 = std::get<0>(i0x);
        std::vector<float> i012 = std::get<1>(i0x);
    } else if (SHADING_MODEL == PHONG) {
        std::tuple nx = EdgeInterpolate(p0.y, normal0.x, p1.y, normal1.x, p2.y, normal2.x);
        std::tuple ny = EdgeInterpolate(p0.y, normal0.y, p1.y, normal1.y, p2.y, normal2.y);
        std::tuple nz = EdgeInterpolate(p0.y, normal0.z, p1.y, normal1.z, p2.y, normal2.z);
        nx02 = std::get<0>(nx);
        nx012 = std::get<1>(nx);
        ny02 = std::get<0>(ny);
        ny012 = std::get<1>(ny);
        nz02 = std::get<0>(nz);
        nz012 = std::get<1>(nz);
    }

    // Determine which is left and which is right.
    std::vector<int> x_left;
    std::vector<int> x_right;
    std::vector<float> iz_left;
    std::vector<float> iz_right;
    std::vector<float> i_left;
    std::vector<float> i_right;
    std::vector<float> nx_left;
    std::vector<float> nx_right;
    std::vector<float> ny_left;
    std::vector<float> ny_right;
    std::vector<float> nz_left;
    std::vector<float> nz_right;
    int m = (int) std::floor(x02.size() / 2);
    if (x02[m] < x012[m]) {
        x_left = x02;
        x_right = x012;
        iz_left = iz02;
        iz_right = iz012;
        i_left = i02;
        i_right = i012;
        nx_left = nx02;
        nx_right = nx012;
        ny_left = ny02;
        ny_right = ny012;
        nz_left = nz02;
        nz_right = nz012;
    } else {
        x_left = x012;
        x_right = x02;
        iz_left = iz012;
        iz_right = iz02;
        i_left = i012;
        i_right = i02;
        nx_left = nx012;
        nx_right = nx02;
        ny_left = ny012;
        ny_right = ny02;
        nz_left = nz012;
        nz_right = nz02;

        std::swap(uz_left, uz_right);
        std::swap(vz_left, vz_right);
    }

    // Draw horizontal segments.
    for (int y = p0.y; y <= p2.y; y++) {
        int xl = x_left[y - p0.y];
        int xr = x_right[y - p0.y];

        // Interpolate attributes for this scanline.
        float zl = iz_left[y - p0.y];
        float zr = iz_right[y - p0.y];
        std::vector<float> zscan = Interpolate(xl, zl, xr, zr);

        std::vector<float> iscan;
        std::vector<float> nxscan;
        std::vector<float> nyscan;
        std::vector<float> nzscan;
        if (SHADING_MODEL == GOURAUD) {
            iscan = Interpolate(xl, i_left[y - p0.y], xr, i_right[y - p0.y]);
        } else if (SHADING_MODEL == PHONG) {
            nxscan = Interpolate(xl, nx_left[y - p0.y], xr, nx_right[y - p0.y]);
            nyscan = Interpolate(xl, ny_left[y - p0.y], xr, ny_right[y - p0.y]);
            nzscan = Interpolate(xl, nz_left[y - p0.y], xr, nz_right[y - p0.y]);
        }

        std::vector<float> uzscan = Interpolate(xl, uz_left[y - p0.y], xr, uz_right[y - p0.y]);
        std::vector<float> vzscan = Interpolate(xl, vz_left[y - p0.y], xr, vz_right[y - p0.y]);

        for (int x = xl; x <= xr; x++) {
            float inv_z = zscan[x - xl];
            if (UpdateDepthBufferIfCloser(x, y, zscan[x - xl], depth_buffer)) {
                if (SHADING_MODEL == FLAT) {
                    // Just use the per-triangle intensity
                } else if (SHADING_MODEL == GOURAUD) {
                    intensity = iscan[x - xl];
                } else if (SHADING_MODEL == PHONG) {
                    Vertex vertex = UnProjectVertex(x, y, inv_z);
                    Vertex normal = Vertex(nxscan[x - xl], nyscan[x - xl], nzscan[x - xl]);
                    intensity = ComputeIllumination(vertex, normal, camera, lights);
                }

                float u = uzscan[x - xl] / zscan[x - xl];
                float v = vzscan[x - xl] / zscan[x - xl];
                rgb color = texture.getTexel(u, v);

                PutPixel(data, x, y, MultiplyColor(intensity, color));
            }
        }
    }
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

void RenderModel(uint8_t data[WIDTH*HEIGHT][3], float depth_buffer[WIDTH*HEIGHT], Model model, Camera camera, std::vector<Light> lights, Mat4x4 orientation, Texture& texture) {
    std::vector<Point> projected;
    for (int i = 0; i < model.vertices.size(); i++) {
        projected.push_back(ProjectVertex(model.vertices[i]));
    }
    for (int i = 0; i < model.triangles.size(); i++) {
        RenderTriangle(data, depth_buffer, model.triangles[i], model.vertices, projected, camera, lights, orientation, texture);
    }
}

void RenderScene(uint8_t data[WIDTH*HEIGHT][3], float depth_buffer[WIDTH*HEIGHT], Camera camera, std::vector<Instance> instances, std::vector<Light> lights) {
    Mat4x4 cameraMatrix = MultiplyMM4(Transposed(camera.orientation), MakeTranslationMatrix(Multiply(-1, camera.position)));

    for (int i = 0; i < instances.size(); i++) {
        Mat4x4 transform = MultiplyMM4(cameraMatrix, instances[i].transform);
        Model clipped = TransformAndClip(camera.clipping_planes, instances[i].model, instances[i].scale, transform);
        if (clipped.is_valid) {
            RenderModel(data, depth_buffer, clipped, camera, lights, instances[i].orientation, instances[i].texture);
        }
    }
}

int main() {
    uint8_t data[WIDTH * HEIGHT][3];
    Clear(data);

    // Cube model
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
        Triangle({0, 1, 2}, RED,    {Vertex( 0,  0,  1), Vertex( 0,  0,  1), Vertex( 0,  0,  1)}, {Point(0, 0), Point(1, 0), Point(1, 1)}),
        Triangle({0, 2, 3}, RED,    {Vertex( 0,  0,  1), Vertex( 0,  0,  1), Vertex( 0,  0,  1)}, {Point(0, 0), Point(1, 1), Point(0, 1)}),
        Triangle({4, 0, 3}, GREEN,  {Vertex( 1,  0,  0), Vertex( 1,  0,  0), Vertex( 1,  0,  0)}, {Point(0, 0), Point(1, 0), Point(1, 1)}),
        Triangle({4, 3, 7}, GREEN,  {Vertex( 1,  0,  0), Vertex( 1,  0,  0), Vertex( 1,  0,  0)}, {Point(0, 0), Point(1, 1), Point(0, 1)}),
        Triangle({5, 4, 7}, BLUE,   {Vertex( 0,  0, -1), Vertex( 0,  0, -1), Vertex( 0,  0, -1)}, {Point(0, 0), Point(1, 0), Point(1, 1)}),
        Triangle({5, 7, 6}, BLUE,   {Vertex( 0,  0, -1), Vertex( 0,  0, -1), Vertex( 0,  0, -1)}, {Point(0, 0), Point(1, 1), Point(0, 1)}),
        Triangle({1, 5, 6}, YELLOW, {Vertex(-1,  0,  0), Vertex(-1,  0,  0), Vertex(-1,  0,  0)}, {Point(0, 0), Point(1, 0), Point(1, 1)}),
        Triangle({1, 6, 2}, YELLOW, {Vertex(-1,  0,  0), Vertex(-1,  0,  0), Vertex(-1,  0,  0)}, {Point(0, 0), Point(1, 1), Point(0, 1)}),
        Triangle({1, 0, 5}, PURPLE, {Vertex( 0,  1,  0), Vertex( 0,  1,  0), Vertex( 0,  1,  0)}, {Point(0, 0), Point(1, 0), Point(1, 1)}),
        Triangle({5, 0, 4}, PURPLE, {Vertex( 0,  1,  0), Vertex( 0,  1,  0), Vertex( 0,  1,  0)}, {Point(0, 1), Point(1, 1), Point(0, 0)}),
        Triangle({2, 6, 7}, CYAN,   {Vertex( 0, -1,  0), Vertex( 0, -1,  0), Vertex( 0, -1,  0)}, {Point(0, 0), Point(1, 0), Point(1, 1)}),
        Triangle({2, 7, 3}, CYAN,   {Vertex( 0, -1,  0), Vertex( 0, -1,  0), Vertex( 0, -1,  0)}, {Point(0, 0), Point(1, 1), Point(0, 1)})
    };

    Model cube = Model(vertices, triangles, Vertex(0, 0, 0), sqrt(3));

    Texture wood_texture = Texture("crate-texture.bmp");
    std::vector<Instance> instances = {
        // TODO: scale is passed twice
        Instance(cube, BuildTransformMatrix(Vertex(-1.5, 0, 7), Identity4x4, 0.75), Identity4x4, wood_texture, 0.75),
        Instance(cube, BuildTransformMatrix(Vertex(1.25, 2.5, 7.5), MakeOYRotationMatrix(195)), MakeOYRotationMatrix(195), wood_texture),
        Instance(cube, BuildTransformMatrix(Vertex(1.75, 0, 5), MakeOYRotationMatrix(-30)), MakeOYRotationMatrix(-30), wood_texture)
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

    std::vector<Light> lights = {
        Light(AMBIENT, 0.2, Vertex(0, 0, 0)),
        Light(DIRECTIONAL, 0.2, Vertex(-1, 0, 1)),
        Light(POINT, 0.6, Vertex(-3, 2, -10))
    };

    float depth_buffer[WIDTH*HEIGHT];
    for (int i = 0; i < WIDTH*HEIGHT; i++) {
        depth_buffer[i] = 0;
    }

    RenderScene(data, depth_buffer, camera, instances, lights);

    if (std::getenv("OUT") && write_bmp_file("output.bmp", data, WIDTH, HEIGHT)) {
        std::cout << "Image written successfully." << std::endl;
    }

    return 0;
}