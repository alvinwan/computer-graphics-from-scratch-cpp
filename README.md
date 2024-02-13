# "Computer Graphics from Scratch" in C++
C++ implementation of demos from "Computer Graphics from Scratch".
*Additionally* provides implementations for the "Extending the Raytracer" and
"Extending the Rasterizer" sections, which didn't previously have
implementations.

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