#!/bin/bash

# Run this file from the root of the repository, like the following
#
#   bash tools/benchmark.sh
#
# As all of the paths in the C++ source code and the scripts here are relative
# to the repository root.

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
    g++ demos/$name.cpp -o main.out -std=c++20 -Ofast
    time for i in {1..10}; do ./main.out; done
done;