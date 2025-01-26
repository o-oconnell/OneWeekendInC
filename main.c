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

typedef struct ray ray;
struct ray {
    vec3 origin;
    vec3 direction;
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

vec3 at(Arena *arena, ray *r, double t) {

    // return r->origin + t * r->direction;
}

vec3 *mul(Arena *arena, double t, vec3 *vec) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = t * vec->e[0];
    result->e[1] = t * vec->e[1];
    result->e[2] = t * vec->e[2];

    return result;
}

vec3 *vec_div(Arena *arena, vec3 *vec, double t) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = vec->e[0] / t;
    result->e[1] = vec->e[1] / t;
    result->e[2] = vec->e[2] / t;

    return result;
}


vec3 *plus(Arena *arena, vec3 *a, vec3 *b) {
    
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = a->e[0] + b->e[0];
    result->e[1] = a->e[1] + b->e[1];
    result->e[2] = a->e[2] + b->e[2];

    return result;
}


vec3 *newVec3(Arena *arena, double r, double g, double b) {
    vec3 *newVec = (vec3*) pushArray(arena, sizeof(vec3));
    newVec->e[0] = r;
    newVec->e[1] = g;
    newVec->e[2] = b;

    return newVec;
}

vec3 *minus(Arena * arena, vec3 *a, vec3 *b) {
    return newVec3(arena, a->e[0] - b->e[0], a->e[1] - b->e[1], a->e[2] - b->e[2]);
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


char *vec2str(Arena *arena, vec3 *vec) {
    char *result = pushArray(arena, 100*sizeof(char));
    sprintf(result, "%lf %lf %lf", vec->e[0], vec->e[1], vec->e[2]);

    return result;
}

int main() {
    size_t arena_size = (1 << 30); // 1 GB
    Arena *arena = arenaAlloc(arena_size);

    double aspect_ratio = 16.0/9.0; // width / height
    int image_width = 400;

    int image_height = (int)(image_width / aspect_ratio); // multiply by 9/16ths
    // image_height = (image_height < 1) ? 1 : image_height;

    double focal_length = 1.0;
    double viewport_height = 2.0;
    double viewport_width = viewport_height * (((double)image_width) / image_height);
    fprintf(stderr, "Viewport height is %lf\n", viewport_height);
    fprintf(stderr, "Viewport width is %lf\n", viewport_width);

    vec3 *camera_center = newVec3(arena, 0, 0, 0);

    // vector from left to right edge of VIEWPORT
    vec3 *viewport_u = newVec3(arena, viewport_width, 0, 0);

    // vector from top to bottom of VIEWPORT
    vec3 *viewport_v = newVec3(arena, 0, -viewport_height, 0);
    
    // the amount of delta we go in the VIEWPORT to get one pixel in the ACTUAL IMAGE
    vec3 *pixel_delta_u = vec_div(arena, viewport_u, image_width);
    vec3 *pixel_delta_v = vec_div(arena, viewport_v, image_height);
    
    fprintf(stderr, "Pixel delta u is %lf\n", x(pixel_delta_u));
    fprintf(stderr, "Pixel delta v is %lf\n", y(pixel_delta_v));

// need to do vector minus NEXT

    // camera_center - (0, 0, focal_length) - viewport_u/2 - viewport_v/2

    // everything is relative to the camera center which is at 0,0,0
    // x axis goes right of camera center
    // y axis up
    // negative z axis in the viewing direction. So when z = -1 we are at the viewport
    // (camera is 1 unit away from the viewport in the z direction) 
    // this is called the FOCAL LENGTH!


    // Camera center + (-focal_length) takes us to the center of the viewport
    // Then minus viewport_u / 2 takes us to the left edge of the viewport
    // 

    fprintf(stderr, "viewport_u/2 IS %s\n", vec2str(arena, minus(arena, camera_center, vec_div(arena, viewport_u, 2.0))));



    // vec3* viewport_upper_left = minus(arena, camera_center,
    //     minus(arena, newVec3(arena, 0, 0, focal_length), 
    //         minus(arena, vec_div(arena, viewport_u, 2.0), vec_div(arena, viewport_v, 2.0))
    //     ));
    
    // WRONG:
    // vec3* viewport_upper_left = minus(arena, camera_center, minus(arena, vec_div(arena, viewport_v, 2.0), vec_div(arena, viewport_u, 2.0)));

    // ah ok. viewport_v - viewport_u gives a negative number. Subtracting that negative makes a positive. So x becomes positive. 
    // not what we want
    // YEAH. THATS BECAUSE SUBTRACTION AINT COMMUTATIVE BOY
    // WE SUBTRACT LEFT TO RIGHT

    //  viewport_upper_left = camera_center - vec3(0, 0, focal_length) - viewport_u/2 - viewport_v/2;
    // first we go negative focal length in the z direction
    // then we go negative viewport width / 2 in the x direction
    // then we go DOUBLE negative (positive) viewport width / 2 in the y direction

    // (in terms of PHYSICAL coordinates, this is correct). We have different coordinates for the actual pixels.


    // Should be: negative in the x direction and z direction, positive in the y direction
    vec3 *viewport_upper_left = minus(arena, 
                minus(arena, 
                    minus(arena, 
                        camera_center, 
                        newVec3(arena, 0, 0, focal_length)), 
                    vec_div(arena, viewport_u, 2)), 
                vec_div(arena, viewport_v, 2));

    vec3 *pixel00_loc = plus(arena, viewport_upper_left, mul(arena, 0.5, plus(arena, pixel_delta_u, pixel_delta_v)));

    // fprintf(stderr, "Camera center is %s\n", vec2str(arena, camera_center));
    // fprintf(stderr, "Viewport v is %s\n", vec2str(arena, viewport_v));
    // fprintf(stderr, "Viewport u is %s\n", vec2str(arena, viewport_u));
    fprintf(stderr, "Viewport upper left is %s\n", vec2str(arena, viewport_upper_left));
    

    // pixel 0,0 is inside of the viewport by pixel_delta_u in the x direction and pixel-delta_v in the y direction
    // so we have a border of size pixel_delta_v on top and bottom and pixel_delta_u on left and right
    fprintf(stderr, "Pixel 0,0 location is %s\n", vec2str(arena, pixel00_loc));    


    // camera_center - 
    
    // newVec3(0, 0, focal_length) -  // go out in the Z direction by the focal length to hit the viewport
    
    // vec_div(arena, viewport_u, 2.0) // since we're in the center go left by half of the viewport width
    
    // - vec_div(arena, viewport_v, 2.0); // since we're in the center, go up by half of the viewport length. Subtract since v is negative.

    // // // Image
    // int image_width = 256;
    // int image_height = 256;

    // Render
    printf("P3\n%d %d\n255\n", image_width, image_height);

    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            double r = (double)(i) / (image_width-1);
            double g = (double)(j) / (image_height-1);
            double b = 0.0;

            vec3 *pixel_color = newVec3(arena, r, g, b);

            write_color(pixel_color);
        }
    }
}