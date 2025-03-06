# One Weekend in C

[Ray Tracing in One Weekend](https://raytracing.github.io/) is a classic course which I saw everywhere when I started learning graphics. I read it twice and learned a lot from it, but found a bunch of drawbacks from the implementation in C++: 
- There's a focus on code design and splitting the code into classes that makes it harder to understand, as well as going off on Cpp-esque tangents like shared pointers.
- There's low-value components like an "interval" class that make the code much more complex than necessary.
- If you follow the design of the code (with significant boilerplate) it can make the code a lot harder to debug and box you into following the book's code conventions, which can lead to missing skills you learn by implementing it from scratch.


Run with `rm -rf test.ppm && gcc main.c -lm && ./a.out > test.ppm`



Final scene from Ray Tracing in One Weekend:

![imagepic3](https://github.com/o-oconnell/RayTracer/blob/main/finalscene.png)


More samples:

![imagepic4](https://github.com/o-oconnell/RayTracer/blob/main/finalscene_detailed.png)

With dielectrics:

![imagepic5](https://github.com/o-oconnell/RayTracer/blob/main/finalscene_dielectric.png)

Dielectrics with more samples:

![imagepic6](https://github.com/o-oconnell/RayTracer/blob/main/finalscene_dielectrics_detailed.png)
