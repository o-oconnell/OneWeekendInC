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

vec3 minus_vararg_nopointer(Arena *arena, vec3 *first, ...) {
    // vec3 *result = (vec3*) pushArray(arena, sizeof(vec3));
    vec3 result;

    va_list args;
    va_start(args, first);

    vec3 *vec = first;
    while (vec != NULL) {
        result.e[0] -= vec->e[0];
        result.e[1] -= vec->e[1];
        result.e[2] -= vec->e[2];

        vec = va_arg(args, vec3*);
    }

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



// 
// |x -c|^2 = r^2
// x = o + du
// 
// Substitute for X
// |o +du - c|^2 = r^2
// (o + du - c) * (o + du - c) = r^2
// o^2 + o * du - oc + du*o + du*du - du * c - c*o -c * du +c^2

// o^2 + c^2 + 2(o * du) - oc + du*du - du*c - c*o - c*du
// Completing the square? Yes, we have (o - c)^2 = o^2 + c^2 -2c*o
// (o - c)^2 + 2(o * du) + du*du - du*c - c*du
// (o - c)^2 + 2(o * du) + du*du -2c * du
// d^2(u*u) - 2c * du + (o - c)^2 + 2(o * du)
// d^2(u*u) - 2c * du + 2o*du + (o - c)^2
// d^2(u*u) + 2(-c * du + o * du) + (o - c) ^2
// d^2(u*u) + 2du * (o - c) + (o - c) ^2
// Now u is a vector so we write this as:
// d^2(u * u) + 2d[ u * (o - c)] + (o - c)^2 - r^2 <-- DONT FORGET THE R^2!

// Now the unknown we have is d (length along the ray)
// a = u*u
// b = 2[ u * (o - c) ]
// c = (o - c)^2 - r^2



// || x - c || ^2 = r^2
// x = o + du




int hit_record_num = 0;

int hit_record_max = 2000;

double hit_sphere(Arena *arena, vec3 *center, double radius, ray *r, double ray_tmin, double ray_tmax, hit_record *rec) {
    
    if (hit_record_num < hit_record_max) {
        // fprintf(stderr, "Origin is %s, center is %s\n", vec2str(arena, &r->origin), vec2str(arena, center));

        fprintf(stderr, "Ray with origin %s and direction %s intersected a sphere with center %s and radius %lf\n", vec2str(arena, &r->origin), vec2str(arena, &r->direction), vec2str(arena, center), radius);
    }

    vec3 *oc = minus(arena, &r->origin, center);

    double a = dot(&r->direction, &r->direction); 
    double b = -2.0 * dot(&r->direction, oc);
    double c = dot(oc, oc) - radius * radius;

    double discriminant = b * b - 4 * a * c;


    // if (y(&r->direction) < 0) {
    //     fprintf(stderr, "Y component < 0, vector is %s\n", vec2str(arena, &r->direction));
    // }

    double root = (-b - sqrt(discriminant)) / (2.0*a);

    if (root <= ray_tmin || root >= ray_tmax || isnan(root)) {
        root = (-b + sqrt(discriminant)) / (2.0 * a);

        if (root <= ray_tmin || root >= ray_tmax || isnan(root)) {
            if (y(&r->direction) < 0) {
                // fprintf(stderr, "Returning FALSE\n");
            }
            return false;
        }
    }
    // if (y(&r->direction) < 0) {
    //     fprintf(stderr, "Returning TRUE\n");
    // }

    rec->t = root;

    vec3* atval = at(arena, r, rec->t); 

    // if (hit_record_num > hit_record_max) {

    //     fprintf(stderr, "The point of intersection is %s, the ray direction is %s\n", vec2str(arena, atval), vec2str(arena, &r->direction));
    // }
    
    hit_record_num++;

    vecCopy(arena, &rec->p, atval);

    // normalizes the normal

    

    vec3 *minusval = minus(arena, &rec->p, center);
    
    vec3 *outward_normal = vec_div(arena, minusval, radius);

    fprintf(stderr, "UHHH??\n");
    vec3 *normal;
    if (dot(&r->direction, outward_normal) > 0.0) {
        // Ray is inside sphere
        fprintf(stderr, "RAY INSIDE SPHERE\n");
        normal = mul(arena, -1.0, outward_normal);
    } else {
        fprintf(stderr, "RAY OUTSIDE SPHERE\n");
        normal = outward_normal;
    }

    vecCopy(arena, &rec->normal, normal);



    return discriminant >= 0;
}

double hit_sphere_alt(Arena *arena, vec3 *center, double radius, ray r, double ray_tmin, double ray_tmax, hit_record *rec) {

    // fprintf(stderr, "Sphere center is %s, radius is %lf\n", vec2str(arena, center), radius);
    // fprintf(stderr, "Intersecting vector direction is %s, origin is %s\n", vec2str(arena, &r->direction), vec2str(arena, &r->origin));

    double a = dot(&r.direction, &r.direction);
    // double b = dot(center, &r.direction) + 2 * dot(&r.direction, &r.origin);
    double b = 2 * (dot(&r.direction, &r.origin) - dot(center, &r.direction));

    // fprintf(stderr, "d*v0 = %lf, c*d = %lf\n", dot(&r->direction, &r->origin), dot(center, &r->direction));
    

    vec3 cv0 = minus_nopointer(center, &r.origin);
    // fprintf(stderr, "ray origin is %s, cv0 is %s\n", vec2str(arena, &r->origin), vec2str(arena, &cv0));
    

    double cv0_squared = dot(&cv0, &cv0);
    double c = cv0_squared - radius * radius;

    // fprintf(stderr, "a = %lf, b = %lf, c = %lf\n", a, b, c);


    double root = (-1 * b  - sqrt(b*b - 4 * a * c)) / (2 * a);

    // changing this to 0.01 instead of like 1 unfucked it for some reason
    if (isnan(root) || root < 0.001 || root > ray_tmax) {
        root = (-1 * b + sqrt(b*b - 4 * a * c)) / (2 * a);

        if (isnan(root) || root < 0.001|| root > ray_tmax) {
            // fprintf(stderr, "Returning false\n");
            return false;
        }
    }

    // fprintf(stderr, "ROOT IS %lf\n", root);

    rec->t = root;
    // vec3 *atval = at(arena, r, rec->t);
    vec3 atval = at_nopointer(&r, rec->t);

    // fprintf(stderr, "REC P IS %s\n", vec2str(arena, atval));
    hit_record_num++;

    vecCopy(arena, &rec->p, &atval);

    // vec3 minusval = minus_nopointer(&rec->p, center);
    // // fprintf(stderr, "Radius is %lf\n", radius);
    // vecCopy(arena, &rec->normal, vec_div(arena, &minusval, radius));


    vec3 minusval = minus_nopointer(&rec->p, center);
    
    vec3 *outward_normal = vec_div(arena, &minusval, radius);

    // fprintf(stderr, "UHHH??\n");
    vec3 *normal;
    if (dot(&r.direction, outward_normal) > 0.0) {
        // Ray is inside sphere
        // fprintf(stderr, "RAY INSIDE SPHERE\n");
        normal = mul(arena, -1.0, outward_normal);
        //     vecCopy(arena, &rec->normal, normal);
        // return false;

    } else {
        // fprintf(stderr, "RAY OUTSIDE SPHERE\n");
        normal = outward_normal;
            // vecCopy(arena, &rec->normal, normal);

    }

    vecCopy(arena, &rec->normal, normal);

    // fprintf(stderr, "REC P IS %s, REC NORMAL IS %s\n", vec2str(arena, &rec->p), vec2str(arena, &rec->normal));

    // fprintf(stderr, "Returning true\n");
    
    return true;
}

#define MAX_DIFFUSE_DEPTH 10
#define HITTABLES_LENGTH 4
#define SAMPLES_PER_PIXEL 10

hittable hittables[HITTABLES_LENGTH];

int worldHit(Arena *arena, ray r, double ray_tmin, double ray_tmax, hit_record *rec) {
    hit_record temp_rec;
    int hit_anything = false;
    double closest_so_far = ray_tmax;

    for (int i = 0; i < HITTABLES_LENGTH; ++i) {
        hittable object = hittables[i];
        // fprintf(stderr, "ITERATING ON HITTABLES %d\n", i);

        if (object.type == SPHERE) {
            if (hit_sphere_alt(arena, &object.center, object.radius, r, ray_tmin, ray_tmax, &temp_rec)) {
                // fprintf(stderr, "HIT SPHERE, THE SPHERE NUMBER IS %d\n", i);
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

    // because normals point out ?
    if (dot(on_unit_sphere, normal) <= 0.0) {
        return on_unit_sphere;
    } else {
        return mul(arena, -1.0, on_unit_sphere);
    }
}

vec3 *reflect(Arena *arena, vec3 *v, vec3 *n) {

    return minus(arena, v,  mul(arena, 2*dot(v, n), n));
}

vec3 *reflect_new(Arena *arena, vec3 *incoming, vec3*n) {

    vec3 *L = mul(arena, -1.0, incoming);

    double NL = 2 * dot(n, L);

    vec3 *comp_a = mul(arena, NL, n);

    return minus(arena, comp_a, L);
}

vec3 *vec_mul(Arena *arena, vec3 *a, vec3 *b) {

    return newVec3(arena, a->e[0]*b->e[0], a->e[1]*b->e[1], a->e[2]*b->e[2]);
}

bool near_zero(vec3 *v) {
    double s = 1e-8;
    return (fabs(v->e[0]) < s) && (fabs(v->e[1]) < s) && (fabs(v->e[2]) < s);
}

vec3 *ray_color(Arena *arena, ray r, int depth) {

    if (depth <= 0) {
        // fprintf(stderr, "------------reached final bounce, returning zero vector\n");
        return newVec3(arena, 0, 0, 0); 
    }

    hit_record rec;
    if (worldHit(arena, r, 0.001, DBL_MAX, &rec)) {
        // in this case since rec.p is always 0,0,0 the direction 
        // is actually same as the ray vector

        // we're facing negative z so I think z doesn't matter?
        // we need |x| >= 1 and/or |y|>=1 to make it out of the sphere though.

        // If we always intersect, that's what causes the black color
        // -0.004444 0.995556 -1.000000
        if (rec.material == LAMBERTIAN) {
            vec3 *random_uv = random_unit_vector(arena);

            // QUESTION: WHY DOES IT MAKE ALL BLACK IF THIS IS A MINUS? 
            // ^TODO: ANSWER THIS QUESTION SO YOU CAN UNDERSTAND YOUR SHIT RETARD
            // ^From plotting a couple points: SEEMS to be the case that minus results in a vector
            // that's always inside the sphere. So always intersects and gives us black.

            vec3 *direction = plus(arena, random_uv, &rec.normal);
            vec3 *minus_direction = plus(arena, random_uv, &rec.normal);


            // fprintf(stderr, "Random unit vector is %s, REC.P is %s, normal is %s, IF I ADD, new direction is %s, IF I SUBTRACT, new direction is %s\n", vec2str(arena, random_uv), vec2str(arena, &rec.p), vec2str(arena, &rec.normal), vec2str(arena, direction), vec2str(arena, minus_direction));

            if (near_zero(direction)) {
                direction = &rec.normal;
            }
            
            double direction_magnitude = sqrt(x(direction) * x(direction) + y(direction) * y(direction) + z(direction) * z(direction));

            ray ray_new = newRay_noArena(rec.p, *direction);

            // fprintf(stderr, "New ray origin %s, direction %s\n", vec2str(arena, &ray_new->origin), vec2str(arena,  &ray_new->direction));

            vec3 *result = vec_mul(arena, &rec.albedo, ray_color(arena, ray_new, depth -1)); 

            // fprintf(stderr, "INTERSECT! Direction is %s, normal is %s, result is %s, DIRECTION MAGNITUDE IS %lf\n", vec2str(arena, direction), vec2str(arena, &rec.normal), vec2str(arena, result), direction_magnitude);
            return result;
        } else {
            vec3 *reflected = reflect_new(arena, &r.direction, &rec.normal);
            
            ray scattered = newRay_noArena(rec.p, *reflected);

            // fprintf(stderr, "Incoming ray origin is %s, direction is %s. Reflected ray origin is %s and direction is %s. Normal ray direction is %s\n", vec2str(arena, &r->origin), vec2str(arena, &r->direction), vec2str(arena, &rec.p), vec2str(arena, reflected), vec2str(arena, &rec.normal));
            // vec3* attenuation_albedo = newVec3(arena, 0.8, 0.8, 0.8);

            if (hit_record_num < hit_record_max) {
                // fprintf(stderr, "Scattered ray origin %s, direction %s\n", vec2str(arena, &scattered->origin), vec2str(arena,  &scattered->direction));
            }

            return vec_mul(arena, &rec.albedo, ray_color(arena, scattered, depth -1));
        }
    }

    vec3 *unit_direction = unit_vector(arena, &r.direction);
    double a = 0.5 * (y(unit_direction) + 1.0);
    
    // makes the blue background. lerp.
    vec3 *result =  plus_vararg(arena, mul(arena, (1.0 - a), newVec3(arena, 1.0, 1.0, 1.0)), mul(arena, a, newVec3(arena, 0.5, 0.7, 1.0)), NULL);

    // fprintf(stderr, "NO INTERSECT! Result is %s\n", vec2str(arena, result));

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

int main() {
    srand(time(NULL));

    size_t arena_size = (1 << 30); // 1 GB
    Arena *arena = arenaAlloc(arena_size);

    vec3 *center2 = newVec3(arena, 0, -100.5, -1);
    vecCopy(arena, &hittables[0].center, center2);
    hittables[0].radius = 100;
    hittables[0].type = SPHERE;
    hittables[0].material = LAMBERTIAN;
    vec3 *albedo0 = newVec3(arena, 0.8, 0.8, 0.0);
    vecCopy(arena, &hittables[0].albedo, albedo0);

    vec3 *center = newVec3(arena, 0, 0, -1.2);
    vecCopy(arena, &hittables[1].center, center);
    hittables[1].radius = 0.5;
    hittables[1].type = SPHERE;
    hittables[1].material = LAMBERTIAN;

    vec3 *albedo = newVec3(arena, 0.1, 0.2, 0.5);
    vecCopy(arena, &hittables[1].albedo, albedo);

    // // // // left
    vec3 *center3 = newVec3(arena, -1.0, 0, -1);
    vecCopy(arena, &hittables[2].center, center3);
    hittables[2].radius = 0.5;
    hittables[2].type = SPHERE;
    hittables[2].material = METAL;

    vec3 *albedo3 = newVec3(arena, 0.8, 0.8, 0.8);
    vecCopy(arena, &hittables[2].albedo, albedo3);

    // // // // right 
    vec3 *center4 = newVec3(arena, 1.0, 0, -1);
    vecCopy(arena, &hittables[3].center, center4);
    hittables[3].radius = 0.5;
    hittables[3].type = SPHERE;
    hittables[3].material = METAL;

    vec3 *albedo4 = newVec3(arena, 0.8, 0.6, 0.2);
    vecCopy(arena, &hittables[3].albedo, albedo4);


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
    int image_width = 400;

    int image_height = (int)(image_width / aspect_ratio); 

    vec3 *lookfrom = newVec3(arena, 0, 0, 0); // why does negative Z values go further out? 
    vec3 *lookat = newVec3(arena, 0, 0, -1);
    vec3 *vup = newVec3(arena, 0, 1, 0);


    // vec3 *camera_center = newVec3(arena, 0, 0, 0);
    vec3 *camera_center = lookfrom;

    double focal_length = length(arena, minus(arena, lookfrom, lookat));

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

    vec3 *w = unit_vector(arena, minus(arena, lookfrom, lookat));

    // right hand rule with vup and pointing towards the z direction
    // gives a vector that points to the right
    vec3 *u = unit_vector(arena, cross(arena, vup, w));

    // right hand rule with the vector pointing in the z direction and 
    // the vector pointing right gives a vector that is similar to vup 
    // but instead of pointing straight up points up "in the plane" normal to 
    // the camera direction
    vec3 *v = cross(arena, w, u);

    vec3 *viewport_u = mul(arena, viewport_width, u);
    
    // why multiply by -1? 
    vec3 *viewport_v = mul(arena, viewport_height, mul(arena, -1.0, v)); 


    // vec3 *viewport_u = newVec3(arena, viewport_width, 0, 0);
    // vec3 *viewport_v = newVec3(arena, 0, -viewport_height, 0);
    
    vec3 *pixel_delta_u = vec_div(arena, viewport_u, image_width);
    vec3 *pixel_delta_v = vec_div(arena, viewport_v, image_height);
    
    fprintf(stderr, "Pixel delta u is %lf\n", x(pixel_delta_u));
    fprintf(stderr, "Pixel delta v is %lf\n", y(pixel_delta_v));


    fprintf(stderr, "viewport_u/2 IS %s\n", vec2str(arena, minus(arena, camera_center, vec_div(arena, viewport_u, 2.0))));


    // mul(arena, focal_length, w) - multiplies unit vector pointing in the Z direction by the actual length
    // viewport_u / 2 = half of the width
    // viewport_v / 2 = half of the height
    vec3 *viewport_upper_left = minus_vararg(arena, camera_center, mul(arena, focal_length, w), vec_div(arena, viewport_u, 2), vec_div(arena, viewport_v, 2), NULL);

    // vec3 *viewport_upper_left = minus_vararg(arena, camera_center, newVec3(arena, 0, 0, focal_length), vec_div(arena, viewport_u, 2), vec_div(arena, viewport_v, 2), NULL);

    vec3 *arg2 = mul(arena, 0.5, plus(arena, pixel_delta_u, pixel_delta_v));
    vec3 pixel00_loc = plus_nopointer(viewport_upper_left, arg2);

    fprintf(stderr, "Viewport upper left is %s\n", vec2str(arena, viewport_upper_left));
    fprintf(stderr, "Pixel 0,0 location is %s\n", vec2str(arena, &pixel00_loc));
    printf("P3\n%d %d\n255\n", image_width, image_height);

    int samples_per_pixel = 10;
    double scale = 1.0/samples_per_pixel;

    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            vec3 pixel_color = newVec3_noArena(0, 0, 0);

            int arg_count = 3;
            vec3 arg2 = mul_nopointer(i, pixel_delta_u); 
            vec3 arg3 = mul_nopointer(j, pixel_delta_v);

            vec3 pixel_center = plus_vararg_nopointer(arg_count, pixel00_loc, arg2, arg3);

            if (hit_record_num < hit_record_max) {
                // fprintf(stderr, "Pixel center is %s (j = %d, i = %d)\n", vec2str(arena, pixel_center), j, i);
            }

            for (int sample = 0; sample < samples_per_pixel; sample++) {
                vec3 ray_direction = minus_nopointer(&pixel_center, camera_center);
                ray r = newRay_noArena(*camera_center, ray_direction);

                if (hit_record_num < hit_record_max) {
                    // fprintf(stderr, "IN MAIN Origin is %s\n", vec2str(arena, &r->origin));
                }

                pixel_color = plus_nopointer_onePointerArgs(pixel_color, ray_color(arena, r, MAX_DIFFUSE_DEPTH));
            }

            write_color(mul_nopointer_noPointerArgs(scale, pixel_color));
        }
    }
}