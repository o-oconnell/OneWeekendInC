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
    LAMBERTIAN,
    DIELECTRIC
};

typedef struct hit_record hit_record;
struct hit_record {
    vec3 p;
    vec3 normal;
    vec3 albedo;
    double t;
    double refraction_index;

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
    double refraction_index;

    enum hittable_type type;
    enum hittable_material material;
};

#define MAX_DIFFUSE_DEPTH 50
#define HITTABLES_LENGTH 500

hittable hittables[HITTABLES_LENGTH];

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

ray newRay(vec3 origin, vec3 direction) {
    ray result;
    result.origin = origin;
    result.direction = direction;

    return result;
}

vec3 newVec3(double r, double g, double b) {
    vec3 newVec;

    newVec.e[0] = r;
    newVec.e[1] = g;
    newVec.e[2] = b;

    return newVec;
}

double length_squared(vec3 v) {
    return v.e[0]*v.e[0] + v.e[1]*v.e[1] + v.e[2]*v.e[2];
}

double length(vec3 v) {
    return sqrt(length_squared(v));
}

vec3 mul(double t, vec3 vec) {
    vec3 result;

    result.e[0] = t * vec.e[0];
    result.e[1] = t * vec.e[1];
    result.e[2] = t * vec.e[2];

    return result;
}

vec3 vec_div(vec3 vec, double t) {
    vec3 result;

    result.e[0] = vec.e[0] / t;
    result.e[1] = vec.e[1] / t;
    result.e[2] = vec.e[2] / t;

    return result;
}

vec3 unit_vector(vec3 v) {
  return vec_div(v, length(v));
}

vec3 plus(vec3 a, vec3 b) { 
    vec3 result;

    result.e[0] = a.e[0] + b.e[0];
    result.e[1] = a.e[1] + b.e[1];
    result.e[2] = a.e[2] + b.e[2];

    return result;
}

vec3 at(ray r, double t) {
    vec3 mulResult = mul(t, r.direction);

    return plus(r.origin, mulResult);
}

vec3 plus_vararg(int count, vec3 first, ...) {
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

    va_end(args);
    return result;
}

vec3 minus_vararg(int count, vec3 first, ...) {
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

    va_end(args);
    return result;
}

vec3 minus(vec3 a, vec3 b) {

    vec3 result;

    result.e[0] = a.e[0] - b.e[0];
    result.e[1] = a.e[1] - b.e[1];
    result.e[2] = a.e[2] - b.e[2];

    return result;
}

double x(vec3 v) {
    return v.e[0];
}

double y(vec3 v) {
    return v.e[1];
}

double z(vec3 v) {
    return v.e[2];
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
    char *result = malloc(100*sizeof(char));
    sprintf(result, "%lf %lf %lf", vec.e[0], vec.e[1], vec.e[2]);

    return result;
}

double dot(vec3 u, vec3 v) {
    return u.e[0] * v.e[0] 
            + u.e[1] * v.e[1]
            + u.e[2] * v.e[2];
}

double hit_sphere_alt(vec3 center, double radius, ray r, double ray_tmin, double ray_tmax, hit_record *rec) {
    double a = dot(r.direction, r.direction);
    double b = 2 * (dot(r.direction, r.origin) - dot(center, r.direction));
    vec3 cv0 = minus(center, r.origin);    

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
    rec->p = at(r, rec->t);

    vec3 outward_normal = vec_div(minus(rec->p, center), radius);
    vec3 normal;
    if (dot(r.direction, outward_normal) > 0.0) {
        rec->front_face = 1;
        normal = mul(-1.0, outward_normal);
    } else {
        rec->front_face = 0;
        normal = outward_normal;
    }

    rec->normal = normal;
    return true;
}

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
                rec->refraction_index = object.refraction_index;
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

vec3 random_vec3() {
    return newVec3(random_double(), random_double(), random_double());
}

vec3 random_range_vec3(double min, double max) {
    return newVec3(random_double_range(min,max), random_double_range(min,max), random_double_range(min,max));
}

vec3 random_unit_vector() {
    while (true) {
        vec3 p = random_range_vec3( -1.0, 1.0);
        double lensq = length_squared(p);
        
        if (1e-160 <= lensq <= 1) {
            return vec_div(p, sqrt(lensq));
        }
    }
}

vec3 reflect(vec3 incoming, vec3 n) {
    vec3 L = mul( -1.0, incoming);
    double NL = 2 * dot(n, L);

    vec3 comp_a = mul(NL, n);
    return minus(comp_a, L);
}

vec3 refract(vec3 uv, vec3 n, double etai_over_etat) {
    double uvn = dot(mul(-1.0, uv), n);
    double cos_theta = uvn;

    if (cos_theta > 1.0) {
        cos_theta = 1.0;
    }

    vec3 r_out_perp = mul(etai_over_etat, plus(uv, mul(cos_theta, n)));

    double r_out_parallel_mul = -1.0 * sqrt(fabs(1.0 - length_squared(r_out_perp)));
    vec3 r_out_parallel = mul(r_out_parallel_mul, n);

    return plus(r_out_perp, r_out_parallel);
}

vec3 vec_mul(vec3 a, vec3 b) {
    return newVec3(a.e[0]*b.e[0], a.e[1]*b.e[1], a.e[2]*b.e[2]);
}

bool near_zero(vec3 *v) {
    double s = 1e-8;
    return (fabs(v->e[0]) < s) && (fabs(v->e[1]) < s) && (fabs(v->e[2]) < s);
}

vec3 ray_color(ray r, int depth) {

    if (depth <= 0) {
        return newVec3(0, 0, 0); 
    }

    hit_record rec;
    if (worldHit(r, 0.001, DBL_MAX, &rec)) {
        if (rec.material == LAMBERTIAN) {
            vec3 random_uv = random_unit_vector();
            vec3 direction = plus(random_uv, rec.normal);
            vec3 minus_direction = plus(random_uv, rec.normal);

            if (near_zero(&direction)) {
                direction = rec.normal;
            }
            
            double direction_magnitude = sqrt(x(direction) * x(direction) + y(direction) * y(direction) + z(direction) * z(direction));

            ray ray_new = newRay(rec.p, direction);
            return vec_mul(rec.albedo, ray_color(ray_new, depth -1)); 
        } else if (rec.material == DIELECTRIC) {
            vec3 attenuation = newVec3(1.0, 1.0, 1.0);
            double ri = rec.front_face ?  rec.refraction_index : (1.0/rec.refraction_index);

            vec3 unit_direction = unit_vector(r.direction);
            double cos_theta = dot(mul(-1.0, unit_direction), rec.normal);
            if (cos_theta > 1.0) cos_theta = 1.0;
            double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

            vec3 direction;
            if (ri * sin_theta > 1.0) {
                direction = reflect(unit_direction, rec.normal);
            } else {
                direction = refract(unit_direction, rec.normal, ri);
            }

            ray scattered = newRay(rec.p, direction);
            return ray_color(scattered, depth - 1);
        } else {
            vec3 reflected = reflect(r.direction, rec.normal);
            ray scattered = newRay(rec.p, reflected);

            return vec_mul(rec.albedo, ray_color(scattered, depth -1));
        }
    }

    vec3 unit_direction = unit_vector(r.direction);
    double a = 0.5 * (y(unit_direction) + 1.0);
    
    return plus(mul((1.0 - a), newVec3(1.0, 1.0, 1.0)), mul(a, newVec3(0.5, 0.7, 1.0)));
}

const double pi = 3.1415926535897932385;

double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

vec3 cross(vec3 u, vec3 v) {
    return newVec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                          u.e[2] * v.e[0] - u.e[0] * v.e[2],
                          u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

int main() {
    srand(time(NULL));

    size_t arena_size = (1 << 30); // 1 GB
    Arena *arena = arenaAlloc(arena_size);

    hittables[0].center = newVec3(0, -1000, 0);
    hittables[0].radius = 1000;
    hittables[0].type = SPHERE;
    hittables[0].material = LAMBERTIAN;
    hittables[0].albedo = newVec3(0.5, 0.5, 0.5);

    int hittables_index = 1;

    for (int a = -11; a < 11; ++a) {
        for (int b = -11; b < 11; ++b) {
            double choose_mat = random_double();

            vec3 center = newVec3(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if (length(minus(center, newVec3(4, 0.2, 0))) > 0.9) {
                hittables[hittables_index].center = center;
                hittables[hittables_index].radius = 0.2;

                if (choose_mat < 0.8) {
                    vec3 albedo = vec_mul(random_vec3(), random_vec3());

                    hittables[hittables_index].material = LAMBERTIAN;
                    hittables[hittables_index].albedo = albedo;
                } else if (choose_mat < 0.95) {
                    vec3 albedo = random_range_vec3(0.5, 1);
                    
                    hittables[hittables_index].material = METAL;
                    hittables[hittables_index].albedo = albedo;
                } else {

                    hittables[hittables_index].refraction_index = 1.5;
                    hittables[hittables_index].material = DIELECTRIC;
                    hittables[hittables_index].albedo = newVec3(0.7, 0.6, 0.5);
                }

                hittables_index++;
            }
        }
    }

    hittables[hittables_index].center = newVec3(0, 1, 0);
    hittables[hittables_index].radius = 1;
    hittables[hittables_index].material = DIELECTRIC;
    hittables[hittables_index].refraction_index = 1.5;
    hittables[hittables_index].albedo = newVec3(0.7, 0.6, 0.5);

    hittables_index++;

    hittables[hittables_index].center = newVec3(-4, 1, 0);
    hittables[hittables_index].radius = 1;
    hittables[hittables_index].material = LAMBERTIAN;
    hittables[hittables_index].albedo = newVec3(0.4, 0.2, 0.1);

    hittables_index++;

    hittables[hittables_index].center = newVec3(4, 1, 0);
    hittables[hittables_index].radius = 1;
    hittables[hittables_index].material = METAL;
    hittables[hittables_index].albedo = newVec3(0.7, 0.6, 0.5);
    
    double aspect_ratio = 16.0/9.0; 
    int image_width = 1200;
    int image_height = (int)(image_width / aspect_ratio); 

    vec3 lookfrom = newVec3(13,2,3);
    vec3 lookat = newVec3(0, 0, 0);
    vec3 vup = newVec3(0, 1, 0);
    vec3 camera_center = lookfrom;

    double focal_length = length(minus(lookfrom, lookat));
    double vfov = 90;
    double theta = degrees_to_radians(vfov);
    double h = tan(theta/2);
    double viewport_height = 2 * h * focal_length;
    double viewport_width = viewport_height * (((double)image_width) / image_height);
    
    fprintf(stderr, "Rendering...\n");
    fprintf(stderr, "Viewport height is %lf\n", viewport_height);
    fprintf(stderr, "Viewport width is %lf\n", viewport_width);

    vec3 w = unit_vector(minus(lookfrom, lookat));
    vec3 u = unit_vector(cross(vup, w));
    vec3 v = cross(w, u);
    vec3 viewport_u = mul(viewport_width, u);
    vec3 viewport_v = mul(viewport_height, mul(-1.0, v)); 
    vec3 pixel_delta_u = vec_div(viewport_u, image_width);
    vec3 pixel_delta_v = vec_div(viewport_v, image_height);
    vec3 viewport_upper_left = minus_vararg(4, camera_center, mul(focal_length, w), vec_div(viewport_u, 2), vec_div(viewport_v, 2));
    vec3 pixel00_loc = plus(viewport_upper_left, mul(0.5, plus(pixel_delta_u, pixel_delta_v)));

    // fprintf(stderr, "Viewport upper left is %s\n", vec2str(arena, &viewport_upper_left));
    // fprintf(stderr, "Pixel 0,0 location is %s\n", vec2str(arena, &pixel00_loc));

    printf("P3\n%d %d\n255\n", image_width, image_height);

    int samples_per_pixel = 500;
    double scale = 1.0/samples_per_pixel;

    for (int j = 0; j < image_height; j++) {
        for (int i = 0; i < image_width; i++) {
            vec3 pixel_color = newVec3(0, 0, 0);

            int arg_count = 3;
            vec3 pixel_center = plus_vararg(arg_count, pixel00_loc, mul(i, pixel_delta_u), mul(j, pixel_delta_v));

            for (int sample = 0; sample < samples_per_pixel; sample++) {
                vec3 ray_direction = minus(pixel_center, camera_center);
                ray r = newRay(camera_center, ray_direction);

                pixel_color = plus(pixel_color, ray_color(r, MAX_DIFFUSE_DEPTH));
            }

            write_color(mul(scale, pixel_color));
        }
    }
}