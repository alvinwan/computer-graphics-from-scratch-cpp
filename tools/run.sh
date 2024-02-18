#!/bin/bash

# Run this file from the root of the repository, like the following
#
#   bash tools/run.sh
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
    "raster-01-lines"
    "raster-02-lines"
    "raster-03-triangle"
    "raster-04-shading"
    "raster-05-perspective"
    "raster-06-cube"
    "raster-07-instances"
    "raster-08-homogenous"
    "raster-09-clipping"
    "raster-10-surface"
    "raster-11-shading"
    "raster-12-texture"
)

# download and convert texture used for raster-12-texture if not already done so
if [[ ! -f crate-texture.bmp ]]; then
   python tools/download_texture.py  # Before running this script, install deps
fi

for name in "${files[@]}"; do
    g++ demos/$name.cpp -o main.out -std=c++20 -Ofast && OUT=true ./main.out && ~/imgcat output.bmp;
done;