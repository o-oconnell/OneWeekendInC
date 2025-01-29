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
#include <stdbool.h>
#include <time.h>
#include <float.h>

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

typedef struct hit_record hit_record;
struct hit_record {
    vec3 p;
    vec3 normal;
    double t;

    int front_face;
};

enum hittable_type {
    SPHERE,
};

typedef struct hittable hittable;
struct hittable {
    vec3 center;
    double radius;

    enum hittable_type type;
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
    dest->e[0] = src->e[0];
    dest->e[1] = src->e[1];
    dest->e[2] = src->e[2];

}

ray* newRay(Arena *arena, vec3 *origin, vec3 *direction) {
    ray *result = (ray*) pushArray(arena, sizeof(ray));

    vecCopy(arena, &result->origin, origin);
    vecCopy(arena, &result->direction, direction);

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
    return plus(arena, &r->origin, mul(arena, t, &r->direction));
}

vec3 *plus_vararg(Arena *arena, vec3 *first, ...) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    va_list args;
    va_start(args, first);

    vec3 *vec = first;
    while (vec != NULL) {
        result->e[0] += vec->e[0];
        result->e[1] += vec->e[1];
        result->e[2] += vec->e[2];

        vec = va_arg(args, vec3*);
    }

    va_end(args);
    return result;
}



vec3 *minus_vararg(Arena *arena, vec3 *first, ...) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    va_list args;
    va_start(args, first);

    vec3 *vec = first;
    while (vec != NULL) {
        result->e[0] -= vec->e[0];
        result->e[1] -= vec->e[1];
        result->e[2] -= vec->e[2];

        vec = va_arg(args, vec3*);
    }

    va_end(args);
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

    int rbyte = (int)(255.999 * r);
    int gbyte = (int)(255.999 * g);
    int bbyte = (int)(255.999 * b);

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

void set_face_normal(Arena *arena, hit_record *rec, ray *r, vec3 *outward_normal) {

    rec->front_face = dot(&r->direction, outward_normal);



    if (rec->front_face) {
        vecCopy(arena, &rec->normal, outward_normal);
    } else {
        vecCopy(arena, &rec->normal, mul(arena, -1.0, outward_normal));
    }
}


double hit_sphere(Arena *arena, vec3 *center, double radius, ray *r, double ray_tmin, double ray_tmax, hit_record *rec) {
    vec3 *oc = minus_vararg(arena, center, &r->origin, NULL);

    double a = dot(&r->direction, &r->direction); 
    double b = -2.0 * dot(&r->direction, oc);
    double c = dot(oc, oc) - radius * radius;

    double discriminant = b * b - 4 * a * c;


    double root = (-b - sqrt(discriminant)) / (2.0*a);

    if (root <= ray_tmin || root >= ray_tmax || isnan(root)) {
        root = (-b + sqrt(discriminant)) / (2.0 * a);

        if (root <= ray_tmin || root >= ray_tmax || isnan(root)) {
            return false;
        }
    }

    rec->t = root;

    vec3* atval = at(arena, r, rec->t); 
    vecCopy(arena, &rec->p, atval);

    vec3 *minusval = minus_vararg(arena, &rec->p, center, NULL);
    vecCopy(arena, &rec->normal, vec_div(arena, minusval, radius));

    return discriminant >= 0;
}


#define MAX_DIFFUSE_DEPTH 10
#define HITTABLES_LENGTH 2
#define SAMPLES_PER_PIXEL 10

hittable hittables[HITTABLES_LENGTH];

int worldHit(Arena *arena, ray *r, double ray_tmin, double ray_tmax, hit_record *rec) {
    hit_record temp_rec;
    int hit_anything = false;
    double closest_so_far = ray_tmax;

    for (int i = 0; i < HITTABLES_LENGTH; ++i) {
        hittable object = hittables[i];

        if (object.type == SPHERE) {
            if (hit_sphere(arena, &object.center, object.radius, r, ray_tmin, ray_tmax, &temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                *rec = temp_rec;
            }
        } else {
            fprintf(stderr, "UNRECOGNIZED OBJECT TYPE\n");
        }
    }

    return hit_anything;
}

double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

double random_double_range(double min, double max) {
    return min + (max - min)*random_double();
}

vec3 *random_vec3(Arena *arena) {
    return newVec3(arena, random_double(), random_double(), random_double());
} 

vec3 *random_range_vec3(Arena *arena, double min, double max) {
    return newVec3(arena, random_double_range(min,max), random_double_range(min,max), random_double_range(min,max));
}

vec3 *random_unit_vector(Arena *arena) {
    while (true) {
        vec3* p = random_range_vec3(arena, -1.0, 1.0);
        double lensq = length_squared(arena, p);
        
        if (1e-160 <= lensq <= 1) {
            return vec_div(arena, p, sqrt(lensq));
        }
    }
}

vec3 *random_on_hemisphere(Arena *arena, vec3 *normal) {
    vec3 *on_unit_sphere = random_unit_vector(arena);

    if (dot(on_unit_sphere, normal) <= 0.0) {
        return on_unit_sphere;
    } else {
        return mul(arena, -1.0, on_unit_sphere);
    }
}

vec3 *reflect(Arena *arena, vec3 *v, vec3 *n) {

    return minus(arena, v,  mul(arena, 2*dot(v, n), n));
}

vec3 *vec_mul(Arena *arena, vec3 *a, vec3 *b) {

    return newVec3(arena, a->e[0]*b->e[0], a->e[1]*b->e[1], a->e[2]*b->e[2]);
}

vec3 *ray_color(Arena *arena, ray* r, int depth) {

    if (depth <= 0) {
        return newVec3(arena, 0, 0, 0); 
    }

    hit_record rec;
    if (worldHit(arena, r, 0, DBL_MAX, &rec)) {

        vec3 *reflected = reflect(arena, &r->direction, &rec.normal);
        ray *scattered = newRay(arena, &rec.p, reflected);
        
      
        vec3* attenuation_albedo = newVec3(arena, 0.8, 0.8, 0.8);

        return vec_mul(arena, attenuation_albedo, ray_color(arena, scattered, depth -1));
    }


    vec3 *unit_direction = unit_vector(arena, &r->direction);
    double a = 0.5 * (y(unit_direction) + 1.0);
    
    vec3 *result =  plus_vararg(arena, mul(arena, (1.0 - a), newVec3(arena, 1.0, 1.0, 1.0)), mul(arena, a, newVec3(arena, 0.5, 0.7, 1.0)), NULL);

    return result;
}

int main() {
    srand(time(NULL));

    size_t arena_size = (1 << 30); // 1 GB
    Arena *arena = arenaAlloc(arena_size);

    vec3 *center2 = newVec3(arena, 0, -100.5, -1);
    vecCopy(arena, &hittables[0].center, center2);
    hittables[0].radius = 100;
    hittables[0].type = SPHERE;


    vec3 *center = newVec3(arena, 0, 0, -1);
    vecCopy(arena, &hittables[1].center, center);
    hittables[1].radius = 0.5;
    hittables[1].type = SPHERE;

    double aspect_ratio = 16.0/9.0; 
    int image_width = 400;

    int image_height = (int)(image_width / aspect_ratio); 

    double focal_length = 1.0;
    double viewport_height = 2.0;
    double viewport_width = viewport_height * (((double)image_width) / image_height);
    fprintf(stderr, "Viewport height is %lf\n", viewport_height);
    fprintf(stderr, "Viewport width is %lf\n", viewport_width);

    vec3 *camera_center = newVec3(arena, 0, 0, 0);

    vec3 *viewport_u = newVec3(arena, viewport_width, 0, 0);
    vec3 *viewport_v = newVec3(arena, 0, -viewport_height, 0);
    
    vec3 *pixel_delta_u = vec_div(arena, viewport_u, image_width);
    vec3 *pixel_delta_v = vec_div(arena, viewport_v, image_height);
    
    fprintf(stderr, "Pixel delta u is %lf\n", x(pixel_delta_u));
    fprintf(stderr, "Pixel delta v is %lf\n", y(pixel_delta_v));


    fprintf(stderr, "viewport_u/2 IS %s\n", vec2str(arena, minus(arena, camera_center, vec_div(arena, viewport_u, 2.0))));


    vec3 *viewport_upper_left = minus_vararg(arena, camera_center, newVec3(arena, 0, 0, focal_length), vec_div(arena, viewport_u, 2), vec_div(arena, viewport_v, 2), NULL);

    vec3 *pixel00_loc = plus_vararg(arena, viewport_upper_left, mul(arena, 0.5, plus(arena, pixel_delta_u, pixel_delta_v)), NULL);

    fprintf(stderr, "Viewport upper left is %s\n", vec2str(arena, viewport_upper_left));
    

    fprintf(stderr, "Pixel 0,0 location is %s\n", vec2str(arena, pixel00_loc));    


    printf("P3\n%d %d\n255\n", image_width, image_height);

    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {


            vec3 *pixel_center = plus_vararg(arena, pixel00_loc, mul(arena, i, pixel_delta_u), mul(arena, j, pixel_delta_v), NULL);
            
            
            vec3 *ray_direction = minus_vararg(arena, pixel_center, camera_center, NULL);


            ray *r = newRay(arena, camera_center, ray_direction);


            vec3 *pixel_color = ray_color(arena, r, MAX_DIFFUSE_DEPTH);

            write_color(pixel_color);
        }
    }
}