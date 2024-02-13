# "Computer Graphics from Scratch" in C++
C++ implementation of demos from "Computer Graphics from Scratch".

Additionally implements features from "Extending the Raytracer" and
"Extending the Rasterizer" sections, which didn't previously have
implementations.

> The C++ demos are written in FP32, as opposed to the original Javascript demos that use FP64 by default. This results in [shadow acne](https://imgur.com/a/ycB69zX) in all the demos after Raycaster 04. To fix this, I added an epsilon reflection ray bias (e.g., start the reflection ray at an epsilon offset along the surface normal). The previews below include this fix.

| Raytracer 01 - Basic | Raytracer 02 - Diffuse | Raytracer 03 - Specular |
|---|---|---|
| ![https://i.imgur.com/ziIToN5.png](https://i.imgur.com/ziIToN5.png) | ![https://i.imgur.com/yaXHtUk.png](https://i.imgur.com/yaXHtUk.png) | ![https://i.imgur.com/bMobaeB.png](https://i.imgur.com/bMobaeB.png) |

| Raytracer 04 - Shadow | Raytracer 05 - Reflection | Raytracer 06 - Camera |
|---|---|---|
| ![https://i.imgur.com/LxoOVys.png](https://i.imgur.com/LxoOVys.png) | ![https://i.imgur.com/RT13thh.png](https://i.imgur.com/RT13thh.png) | ![https://i.imgur.com/B5Tm1CX.png](https://i.imgur.com/B5Tm1CX.png) |

**Bonus demos** - These implement various features from "Extending the Raytracer".

| Raytracer 07 - Triangle | Raytracer 08 - CSG | Raytracer 09 - Optimization |
|---|---|---|
| ![https://i.imgur.com/aieNAAV.png](https://i.imgur.com/aieNAAV.png) | ![https://i.imgur.com/I6g19L4.png](https://i.imgur.com/I6g19L4.png) | ![https://i.imgur.com/cit6dSq.png](https://i.imgur.com/cit6dSq.png) |


## Getting Started

To get started, pick any demo script to build. For example, the first demo for
raycasting is in `raytracer-01.cpp`. To build and run it, use the following

```bash
g++ raytracer-01-basic.cpp -o main.out -std=c++20
./main.out
open output.bmp
```

This gives the following output, which matches the output from the corresponding
[demo on the author's website](https://gabrielgambetta.com/computer-graphics-from-scratch/demos/raytracer-01.html).

![https://i.imgur.com/yEfIoTv.png](https://i.imgur.com/yEfIoTv.png)

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
