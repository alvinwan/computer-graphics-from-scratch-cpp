files=(
    "raytracer-01-basic"
    "raytracer-02-diffuse"
    "raytracer-03-specular"
    "raytracer-04-shadow"
    "raytracer-05-reflection"
    "raytracer-06-camera"
    "raytracer-07-triangles"
    "raytracer-08-csg"
    "raytracer-09-optimization"
)

for name in "${files[@]}"; do
    g++ $name.cpp -o main.out -std=c++20 -Ofast && OUT=true ./main.out && ~/imgcat output.bmp;
done;