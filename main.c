#include <netinet/in.h> 
#include <math.h>
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <sys/socket.h> 
#include <sys/types.h> 
#include <asm-generic/mman.h>
#include <sys/mman.h>
#include <stdarg.h>

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

char *vec2str(Arena *arena, vec3 *vec);

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

void vecCopy(Arena * arena, vec3* dest, vec3 *src) {
    // fprintf(stderr, "DEST IS %d SRC IS %d\n", dest, src);
    // fprintf(stderr, "DEST IS %s, SRC IS %s\n", vec2str(arena, dest), vec2str(arena, src));

    dest->e[0] = src->e[0];
    dest->e[1] = src->e[1];
    dest->e[2] = src->e[2];

    // fprintf(stderr, "DEST IS %s, SRC IS %s\n", vec2str(arena, dest), vec2str(arena, src));

}

ray* newRay(Arena *arena, vec3 *origin, vec3 *direction) {
    ray *result = (ray*) pushArray(arena, sizeof(ray));

    vecCopy(arena, &result->origin, origin);
    vecCopy(arena, &result->direction, direction);

    // fprintf(stderr, "ORIGIN IS %s, DIRECTION IS %s\n", vec2str(arena, origin), vec2str(arena, direction));
    // fprintf(stderr, "NEW RAY ORIGIN IS %s, DIRECTION IS %s\n", vec2str(arena, &result->origin), vec2str(arena, &result->direction));

    return result;
}

vec3 *newVec3(Arena  *arena, double r, double g, double b) {
    vec3 *newVec = (vec3*) pushArray(arena, sizeof(vec3));
    newVec->e[0] = r;
    newVec->e[1] = g;
    newVec->e[2] = b;

    return newVec;
}


double length_squared(Arena *arena, vec3 *v) {
    double result = v->e[0]*v->e[0] + v->e[1]*v->e[1] + v->e[2]*v->e[2];
    // fprintf(stderr, "Length squared result is %lf\n", result);

    return result;  
}

double length(Arena *arena, vec3 *v) {
    return sqrt(length_squared(arena, v));
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

vec3 *unit_vector(Arena *arena, vec3 *v) {
    // vec3 *newVec = (vec3*) pushArray(arena, sizeof(vec3));

    return vec_div(arena, v, length(arena, v));
}


vec3 *plus(Arena *arena, vec3 *a, vec3 *b) {
    
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = a->e[0] + b->e[0];
    result->e[1] = a->e[1] + b->e[1];
    result->e[2] = a->e[2] + b->e[2];

    return result;
}

vec3 *at(Arena *arena, ray *r, double t) {

    // does this assume r->direction is a unit vector ? I think sos
    return plus(arena, &r->origin, mul(arena, t, &r->direction));
    // return r->origin + t * r->direction;
}

// vec3 *plus_vararg(Arena *arena, vec3 *a, ...) {
//     va_list args;
//     va_start(args, a);

//     vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));


//     int result = 0;
//     for (int i = 0; i < count; i++) {
//         result += va_arg(args, int);
//     }

//     va_end(args);
//     return result;
// }

vec3 *plus_vararg(Arena *arena, vec3 *first, ...) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    va_list args;
    va_start(args, first);

    vec3 *vec = first;
    while (vec != NULL) {
        // printf("%s ", vec);
        result->e[0] += vec->e[0];
        result->e[1] += vec->e[1];
        result->e[2] += vec->e[2];

        vec = va_arg(args, vec3*);
    }

    va_end(args);
    // printf("\n");
    return result;
}



vec3 *minus_vararg(Arena *arena, vec3 *first, ...) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    va_list args;
    va_start(args, first);

    vec3 *vec = first;
    while (vec != NULL) {
        // printf("%s ", vec);
        result->e[0] -= vec->e[0];
        result->e[1] -= vec->e[1];
        result->e[2] -= vec->e[2];

        vec = va_arg(args, vec3*);
    }

    va_end(args);
    // printf("\n");
    return result;
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


double dot(vec3 *u, vec3 *v) {

    return u->e[0] * v->e[0] 
            + u->e[1] * v->e[1]
            + u->e[2] * v->e[2];
}


// Book simplifies this to: t^2d⋅d−2td⋅(C−Q)+(C−Q)⋅(C−Q)−r^2=0
// where r is sphere radius, Q is ray origin, d is ray direction, and t is the distance along the ray.
// solve for t.
// it's in the form ax^2 + bc + c = 0
// a = d dot d
// b = 2d dot (C - Q) 
// c = (C - Q) dot (C - Q) - r^2

double hit_sphere(Arena *arena, vec3 *center, double radius, ray *r) {
    vec3 *oc = minus_vararg(arena, center, &r->origin, NULL);

    double a = dot(&r->direction, &r->direction); // length?
    double b = -2.0 * dot(&r->direction, oc);
    double c = dot(oc, oc) - radius * radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-b - sqrt(discriminant)) / (2.0*a);
    }
    // return discriminant >= 0;
}


vec3 *ray_color(Arena *arena, ray* r) {
    double t = hit_sphere(arena, newVec3(arena, 0,0, -1), 0.5, r);
    if (t > 0.0) {

        // This is the normal because it goes from the center of the sphere (0, 0, -1) to the intersection point.
        vec3 *N = unit_vector(arena, minus_vararg(arena, at(arena, r, t), newVec3(arena, 0, 0, -1), NULL));

        // again we add 1 and multiply by 0.5 so the value is in [0, 1] I think
        return mul(arena, 0.5, newVec3(arena, x(N) + 1, y(N) + 1, z(N) + 1));
    }


    // // sphere in the center of the viewport (-1 away from camera)
    // if (hit_sphere(arena, newVec3(arena, 0, 0, -1), 0.5, r)) {
    //     return newVec3(arena, 1, 0, 0);
    // }
    
    // fprintf(stderr, "Ray direction is %s\n", vec2str(arena, &r->direction));

    // this only works because the camera is like 1 unit away from the viewport I think
    // No, not true. It works because we scale it UP to 1.
    vec3 *unit_direction = unit_vector(arena, &r->direction);

    // fprintf(stderr, "unit vector is %s\n", vec2str(arena, unit_direction));
// 
    // Scales the y value from 0 to 1. Reasoning: y component can get as low as -1.0 and as high as 1.0
    // since it's a unit vector. When y components is -1.0, we get 0. When y component is 1.0, we get 1.0;
    double a = 0.5*(y(unit_direction) + 1.0);

    // lerp:
    // when a = 0.0, we get white (1.0, 1.0, 1.0). 
    // when a = 1.1 we get blue (0.5, 0.7, 1.0)

    return plus_vararg(arena, mul(arena, (1.0 - a), newVec3(arena, 1.0, 1.0, 1.0)), mul(arena,  a, newVec3(arena, 0.5, 0.7, 1.0)), NULL);
    // return newVec3(arena, 0,0,0);
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


    // // Should be: negative in the x direction and z direction, positive in the y direction
    // vec3 *viewport_upper_left = minus(arena, 
    //             minus(arena, 
    //                 minus(arena, 
    //                     camera_center, 
    //                     newVec3(arena, 0, 0, focal_length)), 
    //                 vec_div(arena, viewport_u, 2)), 
    //             vec_div(arena, viewport_v, 2));

    vec3 *viewport_upper_left = minus_vararg(arena, camera_center, newVec3(arena, 0, 0, focal_length), vec_div(arena, viewport_u, 2), vec_div(arena, viewport_v, 2), NULL);

    // vec3 *pixel00_loc = plus(arena, viewport_upper_left, mul(arena, 0.5, plus(arena, pixel_delta_u, pixel_delta_v)));

    vec3 *pixel00_loc = plus_vararg(arena, viewport_upper_left, mul(arena, 0.5, plus(arena, pixel_delta_u, pixel_delta_v)), NULL);


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
            // double r = (double)(i) / (image_width-1);
            // double g = (double)(j) / (image_height-1);
            // double b = 0.0;

            // vec3 *pixel_color = newVec3(arena, r, g, b);

            vec3 *pixel_center = plus_vararg(arena, pixel00_loc, mul(arena, i, pixel_delta_u), mul(arena, j, pixel_delta_v), NULL);
            
            // fprintf(stderr, "pixel center is %s\n", vec2str(arena, pixel_center));
            
            vec3 *ray_direction = minus_vararg(arena, pixel_center, camera_center, NULL);

            // fprintf(stderr, "Ray direction HERE is %s\n", vec2str(arena, ray_direction));

            ray *r = newRay(arena, camera_center, ray_direction);

            // printf("After newRay, ray direction is %s\n", vec2str(arena, &r->direction));

            vec3 *pixel_color = ray_color(arena, r);
            // vec3 *pixel_center = plus()

            write_color(pixel_color);
        }
    }
}