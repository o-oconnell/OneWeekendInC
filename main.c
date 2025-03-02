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

enum hittable_material {
    METAL, 
    LAMBERTIAN
};

typedef struct hit_record hit_record;
struct hit_record {
    vec3 p;
    vec3 normal;
    vec3 albedo;
    double t;

    int front_face;
    enum hittable_material material;
    
};

enum hittable_type {
    SPHERE,
};

typedef struct hittable hittable;
struct hittable {
    vec3 center;
    double radius;
    vec3 albedo;

    enum hittable_type type;
    enum hittable_material material;
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

ray newRay_noArena(vec3 origin, vec3 direction) {
    // ray *result = (ray*) pushArray(arena, sizeof(ray));
    ray result;
    result.origin = origin;
    result.direction = direction;

    return result;
}

vec3 *newVec3(Arena  *arena, double r, double g, double b) {
    vec3 *newVec = (vec3*) pushArray(arena, sizeof(vec3));
    newVec->e[0] = r;
    newVec->e[1] = g;
    newVec->e[2] = b;

    return newVec;
}

vec3 newVec3_noArena(double r, double g, double b) {
    vec3 newVec;

    newVec.e[0] = r;
    newVec.e[1] = g;
    newVec.e[2] = b;

    return newVec;
}


double length_squared(Arena *arena, vec3 *v) {
    double result = v->e[0]*v->e[0] + v->e[1]*v->e[1] + v->e[2]*v->e[2];

    return result;  
}

double length_squared_noPointer_noPointerArgs(vec3 v) {
    double result = v.e[0]*v.e[0] + v.e[1]*v.e[1] + v.e[2]*v.e[2];

    return result;  
}

double length(Arena *arena, vec3 *v) {
    return sqrt(length_squared(arena, v));
}

double length_noPointer_noPointerArgs(vec3 v) {
    return sqrt(length_squared_noPointer_noPointerArgs(v));
}

vec3 *mul(Arena *arena, double t, vec3 *vec) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = t * vec->e[0];
    result->e[1] = t * vec->e[1];
    result->e[2] = t * vec->e[2];

    return result;
}

vec3 mul_nopointer(double t, vec3 *vec) {
    vec3 result;

    result.e[0] = t * vec->e[0];
    result.e[1] = t * vec->e[1];
    result.e[2] = t * vec->e[2];

    return result;
}

vec3 mul_nopointer_noPointerArgs(double t, vec3 vec) {
    vec3 result;

    result.e[0] = t * vec.e[0];
    result.e[1] = t * vec.e[1];
    result.e[2] = t * vec.e[2];

    return result;
}

vec3 *vec_div(Arena *arena, vec3 *vec, double t) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = vec->e[0] / t;
    result->e[1] = vec->e[1] / t;
    result->e[2] = vec->e[2] / t;

    return result;
}

vec3 vec_div_noPointer(vec3 *vec, double t) {
    vec3 result;

    result.e[0] = vec->e[0] / t;
    result.e[1] = vec->e[1] / t;
    result.e[2] = vec->e[2] / t;

    return result;
}

vec3 vec_div_noPointer_noPointerArgs(vec3 vec, double t) {
    vec3 result;

    result.e[0] = vec.e[0] / t;
    result.e[1] = vec.e[1] / t;
    result.e[2] = vec.e[2] / t;

    return result;
}

vec3 *unit_vector(Arena *arena, vec3 *v) {
  return vec_div(arena, v, length(arena, v));
}

vec3 unit_vector_noPointer(vec3 v) {
  return vec_div_noPointer_noPointerArgs(v, length_noPointer_noPointerArgs(v));
}


vec3 *plus(Arena *arena, vec3 *a, vec3 *b) {
    
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    result->e[0] = a->e[0] + b->e[0];
    result->e[1] = a->e[1] + b->e[1];
    result->e[2] = a->e[2] + b->e[2];

    return result;
}



vec3 plus_nopointer(vec3 *a, vec3 *b) { 
    vec3 result;

    result.e[0] = a->e[0] + b->e[0];
    result.e[1] = a->e[1] + b->e[1];
    result.e[2] = a->e[2] + b->e[2];

    return result;
}

vec3 plus_nopointer_onePointerArgs(vec3 a, vec3 *b) { 
    vec3 result;

    result.e[0] = a.e[0] + b->e[0];
    result.e[1] = a.e[1] + b->e[1];
    result.e[2] = a.e[2] + b->e[2];

    return result;
}

vec3 plus_nopointer_noPointerArgs(vec3 a, vec3 b) { 
    vec3 result;

    result.e[0] = a.e[0] + b.e[0];
    result.e[1] = a.e[1] + b.e[1];
    result.e[2] = a.e[2] + b.e[2];

    return result;
}


vec3 *at(Arena *arena, ray *r, double t) {
    return plus(arena, &r->origin, mul(arena, t, &r->direction));
}

vec3 at_nopointer(ray *r, double t) {
    vec3 mulResult = mul_nopointer(t, &r->direction);

    return plus_nopointer(&r->origin, &mulResult);
}

vec3 *plus_vararg(Arena *arena, vec3 *first, ...) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));
    result->e[0] = 0;
    result->e[1] = 0;
    result->e[2] = 0;

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

vec3 plus_vararg_nopointer(int count, vec3 first, ...) {
    vec3 result;
    result.e[0] = 0;
    result.e[1] = 0;
    result.e[2] = 0;

    va_list args;
    va_start(args, first);

    vec3 vec = first;

    for (int i = 0; i  < count; ++i) {
        result.e[0] += vec.e[0];
        result.e[1] += vec.e[1];
        result.e[2] += vec.e[2];

        vec = va_arg(args, vec3);
    }
    // while (vec != NULL) {
    //     result.e[0] += vec->e[0];
    //     result.e[1] += vec->e[1];
    //     result.e[2] += vec->e[2];

    //     vec = va_arg(args, vec3);
    // }

    va_end(args);
    return result;
}



vec3 *minus_vararg(Arena *arena, vec3 *first, ...) {
    vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));

    va_list args;
    va_start(args, first);

    // vec3 *vec = first;
    result->e[0] = first->e[0];
    result->e[1] = first->e[1];
    result->e[2] = first->e[2];

    vec3 *vec = va_arg(args, vec3*);

    while (vec != NULL) {
        result->e[0] -= vec->e[0];
        result->e[1] -= vec->e[1];
        result->e[2] -= vec->e[2];

        vec = va_arg(args, vec3*);
    }

    va_end(args);
    return result;
}

vec3 minus_vararg_nopointer(int count, vec3 first, ...) {
    vec3 result;
    result.e[0] = 0;
    result.e[1] = 0;
    result.e[2] = 0;

    va_list args;
    va_start(args, first);

    vec3 vec = first;

    for (int i = 0; i  < count; ++i) {
        result.e[0] -= vec.e[0];
        result.e[1] -= vec.e[1];
        result.e[2] -= vec.e[2];

        vec = va_arg(args, vec3);
    }
    // while (vec != NULL) {
    //     result.e[0] += vec->e[0];
    //     result.e[1] += vec->e[1];
    //     result.e[2] += vec->e[2];

    //     vec = va_arg(args, vec3);
    // }

    va_end(args);
    return result;
}


vec3 *minus(Arena * arena, vec3 *a, vec3 *b) {
    return newVec3(arena, a->e[0] - b->e[0], a->e[1] - b->e[1], a->e[2] - b->e[2]);
}

vec3 minus_nopointer(vec3 *a, vec3 *b) {

    vec3 result;

    result.e[0] = a->e[0] - b->e[0];
    result.e[1] = a->e[1] - b->e[1];
    result.e[2] = a->e[2] - b->e[2];

    return result;
}


vec3 minus_nopointer_noPointerArgs(vec3 a, vec3 b) {

    vec3 result;

    result.e[0] = a.e[0] - b.e[0];
    result.e[1] = a.e[1] - b.e[1];
    result.e[2] = a.e[2] - b.e[2];

    return result;
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

void write_color(vec3 pixel_color) {
    double r = pixel_color.e[0];
    double g = pixel_color.e[1];
    double b = pixel_color.e[2];

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

char *vec2str_noArena(vec3 vec) {
    // char *result = pushArray(arena, 100*sizeof(char));
    char *result = malloc(100*sizeof(char));
    sprintf(result, "%lf %lf %lf", vec.e[0], vec.e[1], vec.e[2]);

    return result;
}


double dot(vec3 u, vec3 v) {

    return u.e[0] * v.e[0] 
            + u.e[1] * v.e[1]
            + u.e[2] * v.e[2];
}


int hit_record_num = 0;

int hit_record_max = 2000;

double hit_sphere_alt(vec3 center, double radius, ray r, double ray_tmin, double ray_tmax, hit_record *rec) {
    double a = dot(r.direction, r.direction);
    double b = 2 * (dot(r.direction, r.origin) - dot(center, r.direction));
    vec3 cv0 = minus_nopointer_noPointerArgs(center, r.origin);    

    double cv0_squared = dot(cv0, cv0);
    double c = cv0_squared - radius * radius;
    double root = (-1 * b  - sqrt(b*b - 4 * a * c)) / (2 * a);

    if (isnan(root) || root < 0.001 || root > ray_tmax) {
        root = (-1 * b + sqrt(b*b - 4 * a * c)) / (2 * a);

        if (isnan(root) || root < 0.001|| root > ray_tmax) {
            return false;
        }
    }

    rec->t = root;
    vec3 atval = at_nopointer(&r, rec->t);

    hit_record_num++;

    rec->p.e[0] = atval.e[0];
    rec->p.e[1] = atval.e[1];
    rec->p.e[2] = atval.e[2];

    vec3 minusval = minus_nopointer_noPointerArgs(rec->p, center);
    vec3 outward_normal = vec_div_noPointer_noPointerArgs(minusval, radius);

    vec3 normal;
    if (dot(r.direction, outward_normal) > 0.0) {
        normal = mul_nopointer_noPointerArgs(-1.0, outward_normal);
    } else {
        normal = outward_normal;
    }

    rec->normal.e[0] = normal.e[0];
    rec->normal.e[1] = normal.e[1];
    rec->normal.e[2] = normal.e[2];
    
    return true;
}

#define MAX_DIFFUSE_DEPTH 10
#define HITTABLES_LENGTH 500
#define SAMPLES_PER_PIXEL 10

hittable hittables[HITTABLES_LENGTH];

int worldHit(ray r, double ray_tmin, double ray_tmax, hit_record *rec) {
    hit_record temp_rec;
    int hit_anything = false;
    double closest_so_far = ray_tmax;

    for (int i = 0; i < HITTABLES_LENGTH; ++i) {
        hittable object = hittables[i];

        if (object.type == SPHERE) {
            if (hit_sphere_alt(object.center, object.radius, r, ray_tmin, closest_so_far, &temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                *rec = temp_rec;

                rec->material = object.material;
                rec->albedo = object.albedo;
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

vec3 random_vec3_noPointer() {
    return newVec3_noArena(random_double(), random_double(), random_double());
} 

vec3 *random_range_vec3(Arena *arena, double min, double max) {
    return newVec3(arena, random_double_range(min,max), random_double_range(min,max), random_double_range(min,max));
}

vec3 random_range_vec3_noPointer(double min, double max) {
    return newVec3_noArena(random_double_range(min,max), random_double_range(min,max), random_double_range(min,max));
}

// vec3 random_range_vec3_noPointer(double min, double max) {
//     return newVec3_noArena(random_double_range(min,max), random_double_range(min,max), random_double_range(min,max));
// }

vec3 random_unit_vector_noArena() {
    while (true) {
        vec3 p = random_range_vec3_noPointer( -1.0, 1.0);
        double lensq = length_squared_noPointer_noPointerArgs(p);
        
        if (1e-160 <= lensq <= 1) {
            return vec_div_noPointer_noPointerArgs(p, sqrt(lensq));
        }
    }
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

vec3 reflect_new(vec3 incoming, vec3 n) {

    vec3 L = mul_nopointer_noPointerArgs( -1.0, incoming);

    double NL = 2 * dot(n, L);

    vec3 comp_a = mul_nopointer_noPointerArgs(NL, n);

    return minus_nopointer_noPointerArgs(comp_a, L);
}

vec3 vec_mul_noPointer(vec3 a, vec3 b) {

    return newVec3_noArena(a.e[0]*b.e[0], a.e[1]*b.e[1], a.e[2]*b.e[2]);
}

bool near_zero(vec3 *v) {
    double s = 1e-8;
    return (fabs(v->e[0]) < s) && (fabs(v->e[1]) < s) && (fabs(v->e[2]) < s);
}

vec3 ray_color(ray r, int depth) {

    if (depth <= 0) {
        return newVec3_noArena(0, 0, 0); 
    }

    hit_record rec;
    if (worldHit(r, 0.001, DBL_MAX, &rec)) {
        if (rec.material == LAMBERTIAN) {
            vec3 random_uv = random_unit_vector_noArena();

            // QUESTION: WHY DOES IT MAKE ALL BLACK IF THIS IS A MINUS? 
            // ^TODO: ANSWER THIS QUESTION SO YOU CAN UNDERSTAND YOUR SHIT RETARD
            // ^From plotting a couple points: SEEMS to be the case that minus results in a vector
            // that's always inside the sphere. So always intersects and gives us black.
            vec3 direction = plus_nopointer_noPointerArgs(random_uv, rec.normal);
            vec3 minus_direction = plus_nopointer_noPointerArgs(random_uv, rec.normal);

            if (near_zero(&direction)) {
                direction = rec.normal;
            }
            
            double direction_magnitude = sqrt(x(&direction) * x(&direction) + y(&direction) * y(&direction) + z(&direction) * z(&direction));

            ray ray_new = newRay_noArena(rec.p, direction);
            vec3 result = vec_mul_noPointer(rec.albedo, ray_color(ray_new, depth -1)); 

            return result;
        } else {
            vec3 reflected = reflect_new(r.direction, rec.normal);
            ray scattered = newRay_noArena(rec.p, reflected);

            return vec_mul_noPointer(rec.albedo, ray_color(scattered, depth -1));
        }
    }

    vec3 unit_direction = unit_vector_noPointer(r.direction);
    double a = 0.5 * (y(&unit_direction) + 1.0);
    
    // makes the blue background. lerp.
    vec3 result =  plus_nopointer_noPointerArgs(mul_nopointer_noPointerArgs((1.0 - a), newVec3_noArena(1.0, 1.0, 1.0)), mul_nopointer_noPointerArgs(a, newVec3_noArena(0.5, 0.7, 1.0)));

    return result;
}

const double pi = 3.1415926535897932385;

double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}


vec3 *cross(Arena *arena, vec3 *u, vec3 *v) {
    return newVec3(arena, u->e[1] * v->e[2] - u->e[2] * v->e[1],
                          u->e[2] * v->e[0] - u->e[0] * v->e[2],
                          u->e[0] * v->e[1] - u->e[1] * v->e[0]);
}

vec3 cross_noPointer(vec3 u, vec3 v) {
    return newVec3_noArena(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                          u.e[2] * v.e[0] - u.e[0] * v.e[2],
                          u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}


int main() {
    srand(time(NULL));

    size_t arena_size = (1 << 30); // 1 GB
    Arena *arena = arenaAlloc(arena_size);

    // vec3 *center2 = newVec3(arena, 0, -100.5, -1);
    // vecCopy(arena, &hittables[0].center, center2);

    // turns black if you set to -1000 - makes sense right? You're inside the sphere. 
    hittables[0].center = newVec3_noArena(0, -1000, 0);
    hittables[0].radius = 1000;
    hittables[0].type = SPHERE;
    hittables[0].material = LAMBERTIAN;
    hittables[0].albedo = newVec3_noArena(0.5, 0.5, 0.5);

    int hittables_index = 1;

    for (int a = -11; a < 11; ++a) {
        for (int b = -11; b < 11; ++b) {
            double choose_mat = random_double();

            vec3 center = newVec3_noArena(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if (length_noPointer_noPointerArgs(minus_nopointer_noPointerArgs(center, newVec3_noArena(4, 0.2, 0))) > 0.9) {
                hittables[hittables_index].center = center;
                hittables[hittables_index].radius = 0.2;

                if (choose_mat < 0.8) {
                    vec3 albedo = vec_mul_noPointer(random_vec3_noPointer(), random_vec3_noPointer());

                    hittables[hittables_index].material = LAMBERTIAN;
                    hittables[hittables_index].albedo = albedo;
                } else {
                    vec3 albedo = random_range_vec3_noPointer(0.5, 1);
                    
                    hittables[hittables_index].material = METAL;
                    hittables[hittables_index].albedo = albedo;
                }

                hittables_index++;
            }
        }
    }


    hittables[hittables_index].center = newVec3_noArena(0, 1, 0);
    hittables[hittables_index].radius = 1;
    hittables[hittables_index].material = METAL;
    hittables[hittables_index].albedo = newVec3_noArena(0.7, 0.6, 0.5);

    hittables_index++;

    hittables[hittables_index].center = newVec3_noArena(-4, 1, 0);
    hittables[hittables_index].radius = 1;
    hittables[hittables_index].material = LAMBERTIAN;
    hittables[hittables_index].albedo = newVec3_noArena(0.4, 0.2, 0.1);

    hittables_index++;

    hittables[hittables_index].center = newVec3_noArena(4, 1, 0);
    hittables[hittables_index].radius = 1;
    hittables[hittables_index].material = METAL;
    hittables[hittables_index].albedo = newVec3_noArena(0.7, 0.6, 0.5);
    



    // vec3 *albedo0 = newVec3(arena, 0.8, 0.8, 0.0);
    // vecCopy(arena, &hittables[0].albedo, albedo0);

    // vec3 *center = newVec3(arena, 0, 0, -1.2);
    // vecCopy(arena, &hittables[1].center, center);

    // hittables[1].center = newVec3_noArena(0, 0, -1.2);
    // hittables[1].radius = 0.5;
    // hittables[1].type = SPHERE;
    // hittables[1].material = LAMBERTIAN;
    // hittables[1].albedo = newVec3_noArena(0.1, 0.2, 0.5);

    // vec3 *albedo = newVec3(arena, 0.1, 0.2, 0.5);
    // vecCopy(arena, &hittables[1].albedo, albedo);

    // // // // // left
    // vec3 *center3 = newVec3(arena, -1.0, 0, -1);
    // vecCopy(arena, &hittables[2].center, center3);

    // hittables[2].center = newVec3_noArena(-1.0, 0, -1);
    // hittables[2].radius = 0.5;
    // hittables[2].type = SPHERE;
    // hittables[2].material = METAL;
    // hittables[2].albedo = newVec3_noArena(0.8, 0.8, 0.8);

    // vec3 *albedo3 = newVec3(arena, 0.8, 0.8, 0.8);
    // vecCopy(arena, &hittables[2].albedo, albedo3);

    // // // // right 
    // vec3 *center4 = newVec3(arena, 1.0, 0, -1);
    // vecCopy(arena, &hittables[3].center, center4);
    // hittables[3].center = newVec3_noArena(1.0, 0, -1);
    // hittables[3].radius = 0.5;
    // hittables[3].type = SPHERE;
    // hittables[3].material = METAL;
    // hittables[3].albedo = newVec3_noArena(0.8, 0.6, 0.2);

    // vec3 *albedo4 = newVec3(arena, 0.8, 0.6, 0.2);
    // vecCopy(arena, &hittables[3].albedo, albedo4);


    // for (int i = 4; i < HITTABLES_LENGTH; ++i) {

    //     vec3 *newCenter = newVec3(arena, i + random_double_range(-1.0, 1.0) - 3, random_double_range(0.0, 1.0), -1);
    //     vecCopy(arena, &hittables[i].center, newCenter);
    //     hittables[i].radius = 0.5;
    //     hittables[i].type = SPHERE;
    //     hittables[i].material = METAL;

    //     vec3 *newAlbedo = newVec3(arena, 0.8, 0.8, 0.8);
    //     vecCopy(arena, &hittables[i].albedo, newAlbedo);

    // }

    double aspect_ratio = 16.0/9.0; 
    int image_width = 1200;

    int image_height = (int)(image_width / aspect_ratio); 

    vec3 lookfrom = newVec3_noArena(13,2,3); // why does negative Z values go further out? 
    vec3 lookat = newVec3_noArena(0, 0, 0);
    vec3 vup = newVec3_noArena(0, 1, 0);


    // vec3 *camera_center = newVec3(arena, 0, 0, 0);
    vec3 camera_center = lookfrom;

    double focal_length = length_noPointer_noPointerArgs(minus_nopointer_noPointerArgs(lookfrom, lookat));

    // double focal_length = 1.0;
    // double viewport_height = 2.0;

    // not sure why this is actually the FOV rather than relating to the height?
    // but it does definitely increase the FOV!
    double vfov = 90;
    double theta = degrees_to_radians(vfov);

    // We have depth of 1 always 
    // so the adjacent has magnitude 1
    // The corner angle of the triangle is theta/2
    // since theta = 90 is the full angle of the camera
    // tan(theta/2) = opposite / adj
    // adj * tan(theta/2) = opposite
    // adj * tan(theta/2) = h
    // Adj = 1, so tan(theta/2) = h
    double h = tan(theta/2);

    // I guess each time we go further inwards in the Z direction the focal length increases,
    // hence multiply by focal length? Otherwise only makes sense to be 2*h.
    double viewport_height = 2 * h * focal_length;

    double viewport_width = viewport_height * (((double)image_width) / image_height);
    fprintf(stderr, "Viewport height is %lf\n", viewport_height);
    fprintf(stderr, "Viewport width is %lf\n", viewport_width);

    vec3 w = unit_vector_noPointer(minus_nopointer_noPointerArgs(lookfrom, lookat));

    // right hand rule with vup and pointing towards the z direction
    // gives a vector that points to the right
    vec3 u = unit_vector_noPointer(cross_noPointer(vup, w));

    // right hand rule with the vector pointing in the z direction and 
    // the vector pointing right gives a vector that is similar to vup 
    // but instead of pointing straight up points up "in the plane" normal to 
    // the camera direction
    vec3 v = cross_noPointer(w, u);

    vec3 viewport_u = mul_nopointer_noPointerArgs(viewport_width, u);
    
    // why multiply by -1? 
    vec3 viewport_v = mul_nopointer_noPointerArgs(viewport_height, mul_nopointer_noPointerArgs(-1.0, v)); 


    // vec3 *viewport_u = newVec3(arena, viewport_width, 0, 0);
    // vec3 *viewport_v = newVec3(arena, 0, -viewport_height, 0);
    
    vec3 pixel_delta_u = vec_div_noPointer_noPointerArgs(viewport_u, image_width);
    vec3 pixel_delta_v = vec_div_noPointer_noPointerArgs(viewport_v, image_height);
    
    fprintf(stderr, "Pixel delta u is %lf\n", x(&pixel_delta_u));
    fprintf(stderr, "Pixel delta v is %lf\n", y(&pixel_delta_v));


    // fprintf(stderr, "viewport_u/2 IS %s\n", vec2str(arena, minus(arena, camera_center, vec_div(arena, viewport_u, 2.0))));


    // mul(arena, focal_length, w) - multiplies unit vector pointing in the Z direction by the actual length
    // viewport_u / 2 = half of the width
    // viewport_v / 2 = half of the height

    vec3 viewport_upper_left = minus_vararg_nopointer(4, camera_center, mul_nopointer_noPointerArgs(focal_length, w), vec_div_noPointer_noPointerArgs(viewport_u, 2), vec_div_noPointer_noPointerArgs(viewport_v, 2));

    // vec3 *viewport_upper_left = minus_vararg(arena, camera_center, newVec3(arena, 0, 0, focal_length), vec_div(arena, viewport_u, 2), vec_div(arena, viewport_v, 2), NULL);

    vec3 arg2 = mul_nopointer_noPointerArgs(0.5, plus_nopointer_noPointerArgs(pixel_delta_u, pixel_delta_v));
    vec3 pixel00_loc = plus_nopointer_noPointerArgs(viewport_upper_left, arg2);
    fprintf(stderr, "Viewport upper left is %s\n", vec2str(arena, &viewport_upper_left));

    fprintf(stderr, "Pixel 0,0 location is %s\n", vec2str(arena, &pixel00_loc));
    printf("P3\n%d %d\n255\n", image_width, image_height);

    int samples_per_pixel = 10;
    double scale = 1.0/samples_per_pixel;

    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            vec3 pixel_color = newVec3_noArena(0, 0, 0);

            int arg_count = 3;
            vec3 arg2 = mul_nopointer_noPointerArgs(i, pixel_delta_u); 
            vec3 arg3 = mul_nopointer_noPointerArgs(j, pixel_delta_v);

            vec3 pixel_center = plus_vararg_nopointer(arg_count, pixel00_loc, arg2, arg3);

            if (hit_record_num < hit_record_max) {
                // fprintf(stderr, "Pixel center is %s (j = %d, i = %d)\n", vec2str(arena, pixel_center), j, i);
            }

            for (int sample = 0; sample < samples_per_pixel; sample++) {
                vec3 ray_direction = minus_nopointer_noPointerArgs(pixel_center, camera_center);
                ray r = newRay_noArena(camera_center, ray_direction);

                if (hit_record_num < hit_record_max) {
                    // fprintf(stderr, "IN MAIN Origin is %s\n", vec2str(arena, &r->origin));
                }

                pixel_color = plus_nopointer_noPointerArgs(pixel_color, ray_color(r, MAX_DIFFUSE_DEPTH));
            }

            write_color(mul_nopointer_noPointerArgs(scale, pixel_color));
        }
    }
}