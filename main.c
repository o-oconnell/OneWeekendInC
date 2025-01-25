// #include <iostream>
#include <netinet/in.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <sys/socket.h> 
#include <sys/types.h> 
#include <asm-generic/mman.h>
#include <sys/mman.h>

typedef struct vec3 vec3;
struct vec3 {
    double e[3];
};

typedef struct Arena Arena;
struct Arena {
    uint64_t pos;
    void* block;
};

Arena *arenaAlloc(size_t size) {
    void *block = mmap(NULL, size, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
    if (block == MAP_FAILED) {
        perror("mmap: ");
        exit(1);
    }

    Arena *arena = (Arena*) block;
    arena->pos = sizeof(Arena);
    arena->block = block;
    return arena;
}

void* pushArray(Arena *arena, size_t size) {
    void *result = (arena->block + arena->pos);
    arena->pos = arena->pos + size;

    return result;
}

vec3 *newVec3(Arena *arena, double r, double g, double b) {
    vec3 *newVec = (vec3*) pushArray(arena, sizeof(vec3));
    newVec->e[0] = r;
    newVec->e[1] = g;
    newVec->e[2] = b;

    return newVec;
}

double x(vec3 *v) {
    return v->e[0];
}

double y(vec3 *v) {
    return v->e[1];
}

double z(vec3 *v) {
    return v->e[2];
}

void write_color(vec3* pixel_color) {
    double r = x(pixel_color);
    double g = y(pixel_color);
    double b = z(pixel_color);

    // Translate the [0,1] component values to the byte range [0,255].
    int rbyte = (int)(255.999 * r);
    int gbyte = (int)(255.999 * g);
    int bbyte = (int)(255.999 * b);

    // Write out the pixel color components.
    printf("%d %d %d\n", rbyte, gbyte, bbyte);
}

int main() {
    size_t arena_size = (1 << 30); // 1 GB
    Arena *arena = arenaAlloc(arena_size);

    // Image
    int image_width = 256;
    int image_height = 256;

    // Render
    printf("P3\n%d %d\n255\n", image_width, image_height);

    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            double r = (double)(i) / (image_width-1);
            double g = (double)(j) / (image_height-1);
            double b = 0.0;

            vec3 *pixel_color = newVec3(arena, r, g, b);

            write_color(pixel_color);

            // int ir = int(255.999 * r);
            // int ig = int(255.999 * g);
            // int ib = int(255.999 * b);

            // printf("%d %d %d\n", ir, ig, ib);

            // double pixel_color = 
        }
    }
}