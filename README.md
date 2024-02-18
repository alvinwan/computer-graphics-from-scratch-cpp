# "Computer Graphics from Scratch" in C++
C++ implementation of demos from "Computer Graphics from Scratch".

Additionally implements features from "Extending the Raytracer" and
"Extending the Rasterizer" sections, which didn't previously have
implementations.

## Raycaster Demos

> The C++ demos are written in FP32, as opposed to the original Javascript demos that use FP64 by default. This results in [shadow acne](https://imgur.com/a/ycB69zX) in all raycasting demos after Raycaster 04. To fix this, I added an epsilon reflection ray bias (e.g., start the reflection ray at an epsilon offset along the surface normal). The previews below include this fix.

| Raytracer 01 - Basic | Raytracer 02 - Diffuse | Raytracer 03 - Specular |
|---|---|---|
| ![https://i.imgur.com/ziIToN5.png](https://i.imgur.com/ziIToN5.png) | ![https://i.imgur.com/yaXHtUk.png](https://i.imgur.com/yaXHtUk.png) | ![https://i.imgur.com/bMobaeB.png](https://i.imgur.com/bMobaeB.png) |

| Raytracer 04 - Shadow | Raytracer 05 - Reflection | Raytracer 06 - Camera |
|---|---|---|
| ![https://i.imgur.com/LxoOVys.png](https://i.imgur.com/LxoOVys.png) | ![https://i.imgur.com/RT13thh.png](https://i.imgur.com/RT13thh.png) | ![https://i.imgur.com/B5Tm1CX.png](https://i.imgur.com/B5Tm1CX.png) |

## Bonus Raytracer Demos

These implement various features from "Extending the Raytracer".

| Raytracer 07 - Triangle | Raytracer 08 - CSG | Raytracer 09 - Optimization |
|---|---|---|
| ![https://i.imgur.com/aieNAAV.png](https://i.imgur.com/aieNAAV.png) | ![https://i.imgur.com/I6g19L4.png](https://i.imgur.com/I6g19L4.png) | ![https://i.imgur.com/cit6dSq.png](https://i.imgur.com/cit6dSq.png) |

## Rasterization Demos

| Raster 01 - Lines | Raster 02 - Lines (Fixed) | Raster 03 - Triangle |
|---|---|---|
| ![https://i.imgur.com/pZXyTU3.png](https://i.imgur.com/pZXyTU3.png) | ![https://i.imgur.com/dALZfnD.png](https://i.imgur.com/dALZfnD.png) | ![https://i.imgur.com/7WuJ6WR.png](https://i.imgur.com/7WuJ6WR.png) |

| Raster 04 - Shading | Raster 05 - Perspective | Raster 06 - Cube |
|---|---|---|
| ![https://i.imgur.com/EWAyk3K.png](https://i.imgur.com/EWAyk3K.png) | ![https://i.imgur.com/GpJoYPp.png](https://i.imgur.com/GpJoYPp.png) | ![https://i.imgur.com/tfZXeZV.png](https://i.imgur.com/tfZXeZV.png) |

| Raster 07 - Instances | Raster 08 - Homogenous | Raster 09 - Clipping |
|---|---|---|
| ![https://i.imgur.com/KrbLL6Z.png](https://i.imgur.com/KrbLL6Z.png) | ![https://i.imgur.com/O0ALrdz.png](https://i.imgur.com/O0ALrdz.png) | ![https://i.imgur.com/qMdhkBA.png](https://i.imgur.com/qMdhkBA.png) |

| Raster 10 - Surface | Raster 11 - Shading | Raster 12 - Texture |
|---|---|---|
| ![https://i.imgur.com/tbswzYt.png](https://i.imgur.com/tbswzYt.png) | ![https://i.imgur.com/Ja4mIuk.png](https://i.imgur.com/Ja4mIuk.png) | ![https://i.imgur.com/G1vMZcN.png](https://i.imgur.com/G1vMZcN.png) |

## Getting Started

To get started, pick any demo script to build. For example, the first demo for
raycasting is in `raytracer-01.cpp`. To build and run it, use the following

```bash
g++ raytracer-01-basic.cpp -o main.out -std=c++20
./main.out
open output.bmp
```

This gives the following output, which matches the output from the original.

| This Repo | Original Reference [on the author's website](https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-01.html) |
|---|---|
| ![https://i.imgur.com/yEfIoTv.png](https://i.imgur.com/yEfIoTv.png) | ![https://i.imgur.com/eJv6tLr.png](https://i.imgur.com/eJv6tLr.png) |

Alternatively, build and run all demos using

```bash
bash run.sh
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
