# "Computer Graphics from Scratch" in C++
Standalone C++ implementation of raytracing and rasterization, without external
libraries or frameworks.
- The demos come from "Computer Graphics from Scratch" ([book](https://gabrielgambetta.com/computer-graphics-from-scratch/)), which provides [Javascript demos](https://github.com/ggambetta/computer-graphics-from-scratch).
- Additionally implements features from "Extending the Raytracer" that weren't implemented in the original source code.

## Raytracer Demos

These demos are all written in FP32, as opposed to the original Javascript demos that use FP64 by default. This results in [shadow acne](https://imgur.com/a/ycB69zX) in all raytracing demos after Raytracer 04. To fix this, I added an epsilon reflection ray bias (e.g., start the reflection ray at an epsilon offset along the surface normal). The previews below include this fix.

| [Raytracer 01 - Basic](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-01-basic.cpp) | [Raytracer 02 - Diffuse](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-02-diffuse.cpp) | [Raytracer 03 - Specular](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-03-specular.cpp) |
|---|---|---|
| ![https://i.imgur.com/ziIToN5.png](https://i.imgur.com/ziIToN5.png) | ![https://i.imgur.com/yaXHtUk.png](https://i.imgur.com/yaXHtUk.png) | ![https://i.imgur.com/bMobaeB.png](https://i.imgur.com/bMobaeB.png) |

| [Raytracer 04 - Shadow](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-04-shadow.cpp) | [Raytracer 05 - Reflection](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-05-reflection.cpp) | [Raytracer 06 - Camera](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-06-camera.cpp) |
|---|---|---|
| ![https://i.imgur.com/LxoOVys.png](https://i.imgur.com/LxoOVys.png) | ![https://i.imgur.com/RT13thh.png](https://i.imgur.com/RT13thh.png) | ![https://i.imgur.com/B5Tm1CX.png](https://i.imgur.com/B5Tm1CX.png) |

## "Extending the Raytracer" Demos

These implement various features from "Extending the Raytracer", which didn't previously have implementations.

| [Raytracer 07 - Triangle](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-07-triangles.cpp) | [Raytracer 08 - CSG](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-08-csg.cpp) | [Raytracer 09 - Optimization](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raytracer-09-optimization.cpp) |
|---|---|---|
| ![https://i.imgur.com/rp5QPvT.png](https://i.imgur.com/rp5QPvT.png) | ![https://i.imgur.com/9kdM0rv.png](https://i.imgur.com/9kdM0rv.png) | ![https://i.imgur.com/wts7Tjg.png](https://i.imgur.com/wts7Tjg.png) |

## Rasterization Demos

The C++ implementations below differ slightly from their Javascript counterparts to improve efficiency.

| [Raster 01 - Lines](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-01-lines.cpp) | [Raster 02 - Lines (Fixed)](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-02-lines.cpp) | [Raster 03 - Triangle](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-03-triangle.cpp) |
|---|---|---|
| ![https://i.imgur.com/pZXyTU3.png](https://i.imgur.com/pZXyTU3.png) | ![https://i.imgur.com/dALZfnD.png](https://i.imgur.com/dALZfnD.png) | ![https://i.imgur.com/7WuJ6WR.png](https://i.imgur.com/7WuJ6WR.png) |

| [Raster 04 - Shading](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-04-shading.cpp) | [Raster 05 - Perspective](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-05-perspective.cpp) | [Raster 06 - Cube](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-06-cube.cpp) |
|---|---|---|
| ![https://i.imgur.com/EWAyk3K.png](https://i.imgur.com/EWAyk3K.png) | ![https://i.imgur.com/GpJoYPp.png](https://i.imgur.com/GpJoYPp.png) | ![https://i.imgur.com/tfZXeZV.png](https://i.imgur.com/tfZXeZV.png) |

| [Raster 07 - Instances](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-07-instances.cpp) | [Raster 08 - Homogenous](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-08-homogenous.cpp) | [Raster 09 - Clipping](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-09-clipping.cpp) |
|---|---|---|
| ![https://i.imgur.com/KrbLL6Z.png](https://i.imgur.com/KrbLL6Z.png) | ![https://i.imgur.com/O0ALrdz.png](https://i.imgur.com/O0ALrdz.png) | ![https://i.imgur.com/qMdhkBA.png](https://i.imgur.com/qMdhkBA.png) |

| [Raster 10 - Surface](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-10-surface.cpp) | [Raster 11 - Shading](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-11-shading.cpp) | [Raster 12 - Texture](https://github.com/alvinwan/computer-graphics-from-scratch-cpp/blob/main/raster-12-texture.cpp) |
|---|---|---|
| ![https://i.imgur.com/tbswzYt.png](https://i.imgur.com/tbswzYt.png) | ![https://i.imgur.com/Ja4mIuk.png](https://i.imgur.com/Ja4mIuk.png) | ![https://i.imgur.com/G1vMZcN.png](https://i.imgur.com/G1vMZcN.png) |

## Getting Started

> **Note**: The instructions below assume you're in the repository root.

To get started, pick any demo script to build. For example, the first demo for
raytracing is in `demos/raytracer-01-basic.cpp`. To build and run it, use the following

```bash
g++ demos/raytracer-01-basic.cpp -o main.out -std=c++20
./main.out
open output.bmp
```

This gives the following output, which matches the output from the original.

| This Repo | Original Reference [on the author's website](https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-01.html) |
|---|---|
| ![https://i.imgur.com/yEfIoTv.png](https://i.imgur.com/yEfIoTv.png) | ![https://i.imgur.com/eJv6tLr.png](https://i.imgur.com/eJv6tLr.png) |

Before running `raster-12-texture.cpp`, download and convert the texture map
using the provided Python script.

```bash
pip install Pillow requests
python tools/download_texture.py
```

Build and run all demos using the following, *after installing Python dependencies above*:

```bash
bash tools/run.sh
```

All scripts in this repository were tested using clang version 15.0.0.

```bash
❯ g++ --version
Apple clang version 15.0.0 (clang-1500.1.0.2.5)
Target: x86_64-apple-darwin23.2.0
Thread model: posix
InstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin
```

<!--
# Develop

The GCC profiler is known to produce empty outputs on a Mac. To run Apple's
profiler, open the `Instruments` app first and use its CLI.

```bash
# open Instruments.app in spotlight
g++ raytracer-01-basic.cpp -o main.out -std=c++20 -Ofast
xctrace record --output 'recording.trace' --template "Time Profiler" --target-stdout - --launch -- ./main.out
```
-->
