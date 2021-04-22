#ifdef MATH_GL
#undef MATH_GL
#define MATH_GL 1
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define PIf 3.14159265358979f
#define HALF_PIf 1.5707963267948966f
#define TAUf 6.283185307179586f
#define SQRT2f 1.4142135623730951f
#define HALF_SQRT2f 0.7071067811865476f
#define SQRT3f 1.7320508075688772f
#define HALF_SQRT3f 0.8660254037844386f
// sqrt(2/3)
#define SQRT_2_3f 0.816496580928f

#include <math.h>

static inline float degrees(float r) {
	return r * (180.0f / PIf);
}
static inline float radians(float r) {
	return r * (PIf / 180.f);
}

// map x from the interval [0, 1] to the interval [a, b]. does NOT clamp.
static inline float lerpf(float x, float a, float b) {
	return x * (b-a) + a;
}

// fractional part of x
static inline float fractf(float x) {
	return x - floorf(x);
}

// opposite of lerp; map x from the interval [a, b] to the interval [0, 1]. does NOT clamp.
static inline float normf(float x, float a, float b) {
	return (x-a) / (b-a);
}

static inline float clampf(float x, float a, float b) {
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

// solve the quadratic equation Ax² + Bx + C = 0, putting the solutions in *x1 and *x2
// if there are 2 solutions, they will be assigned to *x1 and *x2, and true will be returned
// if there is 1 solution, it will be assigned to both *x1 and *x2, and true will be returned (you need to check *x1 == *x2 to detect this case)
// if there are no solutions, *x1 and *x2 will be left intact, and false will be returned
static bool quadratic_equation(float A, float B, float C, float *x1, float *x2) {
	float det = B * B - 4 * A * C;
	if (det < 0) return false;
	float mul = 1.0f / (2*A);
	float sqrt_det = sqrtf(det);
	*x1 = (-B - sqrt_det) * mul;
	*x2 = (-B + sqrt_det) * mul;
	return true;
}

static inline int clampi(int x, int a, int b) {
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static inline i16 clamp_i16(i16 x, i16 a, i16 b) {
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static inline u16 clamp_u16(u16 x, u16 a, u16 b) {
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static inline i32 clamp_i32(i32 x, i32 a, i32 b) {
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static inline u32 clamp_u32(u32 x, u32 a, u32 b) {
	if (x < a) return a;
	if (x > b) return b;
	return x;
}

static inline u8 ndigits_u64(u64 x) {
	u8 ndigits = 1;
	while (x > 9) {
		x /= 10;
		++ndigits;
	}
	return ndigits;
}

// remap x from the interval [from_a, from_b] to the interval [to_a, to_b], NOT clamping if x is outside the "from" interval.
static inline float remapf(float x, float from_a, float from_b, float to_a, float to_b) {
	float pos = (x - from_a) / (from_b - from_a);
	return lerpf(pos, to_a, to_b);
}

static inline float minf(float a, float b) {
	return a < b ? a : b;
}

static inline float maxf(float a, float b) {
	return a > b ? a : b;
}

static inline double maxd(double a, double b) {
	return a > b ? a : b;
}

static inline double mind(double a, double b) {
	return a < b ? a : b;
}

static inline u8 min_u8(u8 a, u8 b) {
	return a < b ? a : b;
}

static inline u8 max_u8(u8 a, u8 b) {
	return a > b ? a : b;
}

static inline u16 min_u16(u16 a, u16 b) {
	return a < b ? a : b;
}

static inline u16 max_u16(u16 a, u16 b) {
	return a > b ? a : b;
}

static inline u32 min_u32(u32 a, u32 b) {
	return a < b ? a : b;
}

static inline u32 max_u32(u32 a, u32 b) {
	return a > b ? a : b;
}

// set *a to the minimum of *a and *b, and *b to the maximum
static inline void sort2_u32(u32 *a, u32 *b) {
	u32 x = *a, y = *b;
	if (x > y) {
		*a = y;
		*b = x;
	}
}

static inline i32 min_i32(i32 a, i32 b) {
	return a < b ? a : b;
}

static inline i32 max_i32(i32 a, i32 b) {
	return a > b ? a : b;
}

static inline u64 min_u64(u64 a, u64 b) {
	return a < b ? a : b;
}

static inline u64 max_u64(u64 a, u64 b) {
	return a > b ? a : b;
}

static inline i64 min_i64(i64 a, i64 b) {
	return a < b ? a : b;
}

static inline i64 max_i64(i64 a, i64 b) {
	return a > b ? a : b;
}

static inline i64 mod_i64(i64 a, i64 b) {
	i64 ret = a % b;
	if (ret < 0) ret += b;
	return ret;
}

static inline i64 abs_i64(i64 x) {
	return x < 0 ? -x : +x;
}

static inline i64 sgn_i64(i64 x) {
	if (x < 0) return -1;
	if (x > 0) return +1;
	return 0;
}

static inline float sgnf(float x) {
	if (x < 0) return -1;
	if (x > 0) return +1;
	return 0;
}

static inline float smoothstepf(float x) {
	if (x <= 0) return 0;
	if (x >= 1) return 1;
	return x * x * (3 - 2 * x);
}

static inline float randf(void) {
	return (float)rand() / (float)((ulong)RAND_MAX + 1);
}

static float rand_gauss(void) {
	// https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
	float U, V;
	do {
		U = randf(), V = randf();
	} while (U == 0 || V == 0);
	return sqrtf(-2 * logf(U)) * cosf(TAUf * V);
}

static u32 rand_u32(void) {
	return ((u32)rand() & 0xfff)
		| ((u32)rand() & 0xfff) << 12
		| ((u32)rand() & 0xff) << 24;
}

static u64 rand_u64(void) {
	return rand_u32() 
		| (u64)rand_u32() << 32;
}

static float rand_uniform(float from, float to) {
	return lerpf(randf(), from, to);
}

#define RAND_SEED_MAX ((u64)0x7FFFFFFFFFFF)
static u64 rand_seed(u64 *seed) {
	// some random numbers
    *seed = (2766354195464186873 * (*seed) + 13309406687124978441u);
    return (*seed>>16) & RAND_SEED_MAX;
}

static float randf_seed(u64 *seed) {
	u64 r = rand_seed(seed);
	return (float)r / (float)(RAND_SEED_MAX + 1);
}

// like randf_seed, but returns in range [-1, +1]
static float rand1_seed(u64 *seed) {
	u64 r = rand_seed(seed);
	float c = (float)RAND_SEED_MAX*.5f;
	return ((float)r - c) / c;
}

static float rand_uniform_seed(u64 *seed, float from, float to) {
	return lerpf(randf_seed(seed), from, to);
}

static float rand_gauss_seed(u64 *seed) {
	float U, V;
	do {
		U = randf_seed(seed), V = randf_seed(seed);
	} while (U == 0 || V == 0);
	return sqrtf(-2 * logf(U)) * cosf(TAUf * V);
}

static float sigmoidf(float x) {
	return 1.0f / (1.0f + expf(-x));
}

// returns ⌈x/y⌉ (x/y rounded up)
static i32 ceildivi32(i32 x, i32 y) {
	if (y < 0) {
		// negating both operands doesn't change the answer
		x = -x;
		y = -y;
	}
	if (x < 0) {
		// truncation is the same as ceiling for negative numbers
		return x / y;
	} else {
		return (x + (y-1)) / y;
	}
}

typedef struct {
	float x, y;
} v2;

static v2 const v2_zero = {0, 0};
static v2 V2(float x, float y) {
	v2 v;
	v.x = x;
	v.y = y;
	return v;
}

static inline v2 v2_add(v2 a, v2 b) {
	return V2(a.x + b.x, a.y + b.y);
}

// a + b * s
static inline v2 v2_add_scaled(v2 a, v2 b, float s) {
	return V2(a.x + b.x * s, a.y + b.y * s);
}

static inline v2 v2_add_const(v2 a, float c) {
	return V2(a.x + c, a.y + c);
}

static inline v2 v2_sub(v2 a, v2 b) {
	return V2(a.x - b.x, a.y - b.y);
}

static inline v2 v2_scale(v2 v, float s) {
	return V2(v.x * s, v.y * s);
}

static inline v2 v2_mul(v2 a, v2 b) {
	return V2(a.x * b.x, a.y * b.y);
}

static inline v2 v2_clamp(v2 x, v2 a, v2 b) {
	return V2(clampf(x.x, a.x, b.x), clampf(x.y, a.y, b.y));
}

static inline float v2_dot(v2 a, v2 b) {
	return a.x * b.x + a.y * b.y;
}

static inline float v2_len(v2 v) {
	return sqrtf(v2_dot(v, v));
}

static inline v2 v2_lerp(float x, v2 a, v2 b) {
	return V2(lerpf(x, a.x, b.x), lerpf(x, a.y, b.y));
}

// rotate v theta radians counterclockwise
static v2 v2_rotate(v2 v, float theta) {
	float c = cosf(theta), s = sinf(theta);
	return V2(
		c * v.x - s * v.y,
		s * v.x + c * v.y
	);
}

static v2 v2_normalize(v2 v) {
	float len = v2_len(v);
	float mul = len == 0.0f ? 1.0f : 1.0f/len;
	return v2_scale(v, mul);
}

static float v2_dist(v2 a, v2 b) {
	return v2_len(v2_sub(a, b));
}

static float v2_dist_squared(v2 a, v2 b) {
	v2 diff = v2_sub(a, b);
	return v2_dot(diff, diff);
}

static void v2_print(v2 v) {
	printf("(%f, %f)\n", v.x, v.y);
}

static v2 v2_rand_unit(void) {
	float theta = rand_uniform(0, TAUf);
	return V2(cosf(theta), sinf(theta));
}

static v2 v2_polar(float r, float theta) {
	return V2(r * cosf(theta), r * sinf(theta));
}

typedef struct {
	float x, y, z;
} v3;

static v3 const v3_zero = {0, 0, 0};

static v3 V3(float x, float y, float z) {
	v3 v;
	v.x = x;
	v.y = y;
	v.z = z;
	return v;
}

static inline v3 v3_from_v2(v2 v) {
	return V3(v.x, v.y, 0);
}

static inline v3 v3_add(v3 a, v3 b) {
	return V3(a.x + b.x, a.y + b.y, a.z + b.z);
}

static inline v3 v3_add_const(v3 a, float c) {
	return V3(a.x + c, a.y + c, a.z + c);
}

// a + b * s
static inline v3 v3_add_scaled(v3 a, v3 b, float s) {
	return V3(a.x + b.x * s, a.y + b.y * s, a.z + b.z * s);
}

// add to only the x component
static inline v3 v3_add_x(v3 a, float dx) {
	return V3(a.x + dx, a.y, a.z);
}
static inline v3 v3_add_y(v3 a, float dy) {
	return V3(a.x, a.y + dy, a.z);
}
static inline v3 v3_add_z(v3 a, float dz) {
	return V3(a.x, a.y, a.z + dz);
}

static inline v3 v3_sub(v3 a, v3 b) {
	return V3(a.x - b.x, a.y - b.y, a.z - b.z);
}

static inline v3 v3_scale(v3 v, float s) {
	return V3(v.x * s, v.y * s, v.z * s);
}


static inline v3 v3_mul(v3 a, v3 b) {
	return V3(a.x * b.x, a.y * b.y, a.z * b.z);
}

static inline v3 v3_lerp(float x, v3 a, v3 b) {
	return V3(lerpf(x, a.x, b.x), lerpf(x, a.y, b.y), lerpf(x, a.z, b.z));
}

static inline float v3_dot(v3 u, v3 v) {
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

static inline v3 v3_cross(v3 u, v3 v) {
	return V3(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
}

static inline float v3_len(v3 v) {
	return sqrtf(v3_dot(v, v));
}

// normalize, then scale
static v3 v3_set_scale(v3 v, float s) {
	float m = s / v3_len(v);
	return v3_scale(v, m);
}

static float v3_dist(v3 a, v3 b) {
	return v3_len(v3_sub(a, b));
}

static inline float v3_dist_squared(v3 a, v3 b) {
	v3 diff = v3_sub(a, b);
	return v3_dot(diff, diff);
}

static v3 v3_normalize(v3 v) {
	float mul = 1.0f / sqrtf(v3_dot(v, v));
	return v3_scale(v, mul);
}

static inline v2 v3_xy(v3 v) {
	return V2(v.x, v.y);
}

static inline void v3_uniform(GLint uniform, v3 v) {
	glUniform3f(uniform, v.x, v.y, v.z);
}

// a point on a unit sphere 
static inline v3 v3_on_sphere(float yaw, float pitch) {
	return V3(cosf(yaw) * cosf(pitch), sinf(pitch), sinf(yaw) * cosf(pitch));
}

static void v3_print(v3 v) {
	printf("(%f, %f, %f)\n", v.x, v.y, v.z);
}

static inline v3 v3_rand(void) {
	return V3(randf(), randf(), randf());
}

static inline v3 v3_rand_seed(u64 *seed) {
	float x = randf_seed(seed);
	float y = randf_seed(seed);
	float z = randf_seed(seed);
	return V3(x, y, z);
}

static inline v3 v3_rand1_seed(u64 *seed) {
	float x = rand1_seed(seed);
	float y = rand1_seed(seed);
	float z = rand1_seed(seed);
	return V3(x, y, z);
}

static v3 v3_rand_unit_seed(u64 *seed) {
	float alpha = acosf(randf_seed(seed) * 2.0f - 1.0f);
	float beta = randf_seed(seed) * TAUf;
	float ca = cosf(alpha), sa = sinf(alpha), cb = cosf(beta), sb = sinf(beta);
	return V3(sa * cb, sa * sb, ca);
}

static v3 v3_rand_unit(void) {
	u64 rng = rand_u64();
	return v3_rand_unit_seed(&rng);
}

static v3 v3_rand_gauss_seed(u64 *seed, float mean, float stddev) {
	float x = rand_gauss_seed(seed) * stddev + mean;
	float y = rand_gauss_seed(seed) * stddev + mean;
	float z = rand_gauss_seed(seed) * stddev + mean;
	return V3(x, y, z);
}

static v3 v3_rand_gauss(float mean, float stddev) {
	u64 seed = rand_u64();
	return v3_rand_gauss_seed(&seed, mean, stddev);
}

// move two vectors towards each other
// x = 0 => don't move them.
// x = 1 => move them both to the center
static void v3_join(v3 *a, v3 *b, float x) {
	v3 center = v3_lerp(0.5f, *a, *b);
	*a = v3_lerp(x, *a, center);
	*b = v3_lerp(x, *b, center);
}

typedef struct {
	float x, y, z, w;
} v4;

static v4 const v4_zero = {0, 0, 0, 0};

static v4 V4(float x, float y, float z, float w) {
	v4 v;
	v.x = x;
	v.y = y;
	v.z = z;
	v.w = w;
	return v;
}

static v4 v4_from_v3(v3 v, float w) {
	return V4(v.x, v.y, v.z, w);
}

static v4 v4_add(v4 a, v4 b) {
	return V4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

static v4 v4_sub(v4 a, v4 b) {
	return V4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

static v4 v4_scale(v4 v, float s) {
	return V4(v.x * s, v.y * s, v.z * s, v.w * s);
}

static v4 v4_scale_xyz(v4 v, float s) {
	return V4(v.x * s, v.y * s, v.z * s, v.w);
}

static v4 v4_lerp(float x, v4 a, v4 b) {
	return V4(lerpf(x, a.x, b.x), lerpf(x, a.y, b.y), lerpf(x, a.z, b.z), lerpf(x, a.w, b.w));
}

static float v4_dot(v4 u, v4 v) {
	return u.x*v.x + u.y*v.y + u.z*v.z + u.w*v.w;
}

// create a new vector by multiplying the respective components of u and v 
static v4 v4_mul(v4 u, v4 v) {
	return V4(u.x * v.x, u.y * v.y, u.z * v.z, u.w * v.w);
}

static float v4_len(v4 v) {
	return sqrtf(v4_dot(v, v));
}

static v4 v4_normalize(v4 v) {
	float len = v4_len(v);
	float mul = len == 0.0f ? 1.0f : 1.0f/len;
	return v4_scale(v, mul);
}

static v3 v4_xyz(v4 v) {
	return V3(v.x, v.y, v.z);
}

static v4 v4_rand(void) {
	return V4(randf(), randf(), randf(), randf());
}

static void v4_uniform(GLint uniform, v4 v) {
	glUniform4f(uniform, v.x, v.y, v.z, v.w);
}

static void v4_print(v4 v) {
	printf("(%f, %f, %f, %f)\n", v.x, v.y, v.z, v.w);
}

typedef struct {
	double x, y;
} v2d;

static v2d V2D(double x, double y) {
	v2d v;
	v.x = x;
	v.y = y;
	return v;
}

// NOTE: matrices are column-major, because that's what they are in OpenGL 

typedef struct {
	float e[4];
} m2;

static m2 const m2_identity = {{
	1, 0,
	0, 1
}};

static void m2_print(m2 m) {
	int i;
	for (i = 0; i < 2; ++i)
		printf("[ %f %f ]\n", m.e[i], m.e[i+2]);
	printf("\n");
}

static m2 M2(float a, float b, float c, float d) {
	m2 ret;
	float *x = ret.e;
	x[0] = a; x[2] = b;
	x[1] = c; x[3] = d;
	return ret;
}

static inline void m2_uniform(GLint u, m2 const *mat) {
	glUniformMatrix2fv(u, 1, GL_FALSE, mat->e);
}

typedef struct {
	float e[9];
} m3;

static m3 m3_identity = {{
	1, 0, 0,
	0, 1, 0,
	0, 0, 1
}};

static void m3_print(m3 m) {
	int i;
	for (i = 0; i < 3; ++i) {
		printf("[ %f %f %f ]\n", m.e[i], m.e[i+3], m.e[i+6]);
	}
	printf("\n");
}

static m3 M3(
	float a, float b, float c,
	float d, float e, float f,
	float g, float h, float i) {
	m3 ret;
	float *x = ret.e;
	x[0] = a; x[3] = b; x[6] = c;
	x[1] = d; x[4] = e; x[7] = f;
	x[2] = g; x[5] = h; x[8] = i;
	return ret;
}

static m3 m3_translate(v2 t) {
	return M3(
		1, 0, t.x,
		0, 1, t.y,
		0, 0, 1
	);
}

static m3 m3_rotate(float theta) {
	float ct = cosf(theta), st = sinf(theta);
	return M3(
		ct, -st, 0,
		st,  ct, 0,
		0,    0, 1
	);
}

static inline void m3_uniform(GLint u, m3 const *mat) {
	glUniformMatrix3fv(u, 1, GL_FALSE, mat->e);
}

static m3 m3_mul(m3 a, m3 b) {
	m3 prod = {0};
	int i, j;
	float *x = prod.e;
	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j, ++x) {
			float *as = &a.e[j];
			float *bs = &b.e[3*i];
			*x = as[0]*bs[0] + as[3]*bs[1] + as[6]*bs[2];
		}
	}
	return prod;
}

typedef struct {
	float e[16];
} m4;

static m4 const m4_identity = {{
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0, 1
}};

static void m4_print(m4 m) {
	int i;
	for (i = 0; i < 4; ++i)
		printf("[ %f %f %f %f ]\n", m.e[i], m.e[i+4], m.e[i+8], m.e[i+12]);
	printf("\n");
}

static m4 M4(
	float a, float b, float c, float d,
	float e, float f, float g, float h,
	float i, float j, float k, float l,
	float m, float n, float o, float p) {
	m4 ret;
	float *x = ret.e;
	x[0] = a; x[4] = b; x[ 8] = c; x[12] = d;
	x[1] = e; x[5] = f; x[ 9] = g; x[13] = h;
	x[2] = i; x[6] = j; x[10] = k; x[14] = l;
	x[3] = m; x[7] = n; x[11] = o; x[15] = p;
	return ret;
}

// see https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations 
static m4 m4_yaw(float yaw) {
	float c = cosf(yaw), s = sinf(yaw);
	return M4(
		c, 0, -s, 0,
		0, 1,  0, 0,
		s, 0,  c, 0,
		0, 0,  0, 1
	);
}

static m4 m4_pitch(float pitch) {
	float c = cosf(pitch), s = sinf(pitch);
	return M4(
		1, 0,  0, 0,
		0, c, -s, 0,
		0, s,  c, 0,
		0, 0,  0, 1
	);
}

// https://en.wikipedia.org/wiki/Translation_(geometry) 
static m4 m4_translate(v3 t) {
	return M4(
		1, 0, 0, t.x,
		0, 1, 0, t.y,
		0, 0, 1, t.z,
		0, 0, 0, 1
	);
}

// multiply m by [v.x, v.y, v.z, 1] 
static v3 m4_mul_v3(m4 m, v3 v) {
	return v3_add(v3_scale(V3(m.e[0], m.e[1], m.e[2]), v.x), v3_add(v3_scale(V3(m.e[4], m.e[5], m.e[6]), v.y),
		v3_add(v3_scale(V3(m.e[8], m.e[9], m.e[10]), v.z), V3(m.e[12], m.e[13], m.e[14]))));
}

/*
4x4 perspective matrix.
fov - field of view in radians, aspect - width:height aspect ratio, z_near/z_far - clipping planes
math stolen from gluPerspective (https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml)
*/
static m4 m4_perspective(float fov, float aspect, float z_near, float z_far) {
	float f = 1.0f / tanf(fov / 2.0f);
	return M4(
		f/aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (z_far+z_near) / (z_near-z_far), (2.0f*z_far*z_near) / (z_near-z_far),
		0, 0, -1, 0
	);
}

// windows.h defines near and far, so let's not use those 
static m4 m4_ortho(float left, float right, float bottom, float top, float near_, float far_) {
	float tx = -(right + left)/(right - left);
	float ty = -(top + bottom)/(top - bottom);
	float tz = -(far_ + near_)/(far_ - near_);
	return M4(
		2.0f / (right - left), 0, 0, tx,
		0, 2.0f / (top - bottom), 0, ty,
		0, 0, -2.0f / (far_ - near_), tz,
		0, 0, 0, 1
	);
}


static m4 m4_mul(m4 a, m4 b) {
	m4 prod = {0};
	int i, j;
	float *x = prod.e;
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j, ++x) {
			float *as = &a.e[j];
			float *bs = &b.e[4*i];
			*x = as[0]*bs[0] + as[4]*bs[1] + as[8]*bs[2] + as[12]*bs[3];
		}
	}
	return prod;
}

static m4 m4_inv(m4 mat) {
	m4 ret;
	float *inv = ret.e;
	float *m = mat.e;

	inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
	inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
	inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
	inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
	inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
	inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
	inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
	inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
	inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
	inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
	inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
	inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
	inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
	inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
	inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
	inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

	float det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0) {
		memset(&ret, 0, sizeof ret);
	} else {
		det = 1 / det;

		for (int i = 0; i < 16; i++)
			inv[i] *= det;
	}

	return ret;
}

static inline void m4_uniform(GLint u, m4 const *mat) {
	glUniformMatrix4fv(u, 1, GL_FALSE, mat->e);
}

typedef struct {
	int x, y;
} v2i;

static v2i V2I(int x, int y) {
	v2i v;
	v.x = x;
	v.y = y;
	return v;
}

static void rgba_u32_to_floats(u32 rgba, float floats[4]) {
	floats[0] = (float)((rgba >> 24) & 0xFF) / 255.f;
	floats[1] = (float)((rgba >> 16) & 0xFF) / 255.f;
	floats[2] = (float)((rgba >>  8) & 0xFF) / 255.f;
	floats[3] = (float)((rgba >>  0) & 0xFF) / 255.f;
}

static u32 floats_to_rgba_u32(float const floats[4]) {
	return (u32)(floats[0] * 255) << 24
	     | (u32)(floats[1] * 255) << 16
	     | (u32)(floats[2] * 255) << 8
	     | (u32)(floats[3] * 255);
}

static v4 rgba_u32_to_v4(u32 rgba) {
	float c[4];
	rgba_u32_to_floats(rgba, c);
	return V4(c[0], c[1], c[2], c[3]);
}

static u32 v3_to_rgba_u32(v3 v) {
	assert(v.x >= 0 && v.x <= 1 && v.y >= 0 && v.y <= 1 && v.z >= 0 && v.z <= 1);
	u32 r = (u32)(v.x * 255);
	u32 g = (u32)(v.y * 255);
	u32 b = (u32)(v.z * 255);
	return r << 24 | g << 16 | b << 8 | 0xFF;
}

static u32 v4_to_rgba_u32(v4 v) {
	assert(v.x >= 0 && v.x <= 1 && v.y >= 0 && v.y <= 1 && v.z >= 0 && v.z <= 1 && v.w >= 0 && v.w <= 1);
	u32 r = (u32)(v.x * 255);
	u32 g = (u32)(v.y * 255);
	u32 b = (u32)(v.z * 255);
	u32 a = (u32)(v.w * 255);
	return r << 24 | g << 16 | b << 8 | a;
}

// returns average of red green and blue components of color
static float rgba_brightness(u32 color) {
	u8 r = (u8)(color >> 24), g = (u8)(color >> 16), b = (u8)(color >> 8);
	return ((float)r+(float)g+(float)b) * (1.0f / 3);
}

static bool rect_contains_point_v2(v2 pos, v2 size, v2 point) {
	float x1 = pos.x, y1 = pos.y, x2 = pos.x + size.x, y2 = pos.y + size.y,
		x = point.x, y = point.y;
	return x >= x1 && x < x2 && y >= y1 && y < y2;
}

static bool centered_rect_contains_point(v2 center, v2 size, v2 point) {
	return rect_contains_point_v2(v2_sub(center, v2_scale(size, 0.5f)), size, point);
}

typedef struct {
	v2 pos, size;
} Rect;

static Rect rect(v2 pos, v2 size) {
	Rect r;
	r.pos = pos;
	r.size = size;
	return r;
}

static Rect rect4(float x1, float y1, float x2, float y2) {
	assert(x2 >= x1);
	assert(y2 >= y1);
	return rect(V2(x1,y1), V2(x2-x1, y2-y1));
}

static Rect rect_centered(v2 center, v2 size) {
	Rect r;
	r.pos = v2_sub(center, v2_scale(size, 0.5f));
	r.size = size;
	return r;
}

static v2 rect_center(Rect r) {
	return v2_add(r.pos, v2_scale(r.size, 0.5f));
}

static bool rect_contains_point(Rect r, v2 point) {
	return rect_contains_point_v2(r.pos, r.size, point);
}

static Rect rect_translate(Rect r, v2 by) {
	return rect(v2_add(r.pos, by), r.size);
}

static float rect_x1(Rect r) { return r.pos.x; }
static float rect_y1(Rect r) { return r.pos.y; }
static float rect_x2(Rect r) { return r.pos.x + r.size.x; }
static float rect_y2(Rect r) { return r.pos.y + r.size.y; }
static float rect_xmid(Rect r) { return r.pos.x + r.size.x * 0.5f; }
static float rect_ymid(Rect r) { return r.pos.y + r.size.y * 0.5f; }

static void rect_coords(Rect r, float *x1, float *y1, float *x2, float *y2) {
	*x1 = r.pos.x;
	*y1 = r.pos.y;
	*x2 = r.pos.x + r.size.x;
	*y2 = r.pos.y + r.size.y;
}

static void rect_print(Rect r) {
	printf("Position: (%f, %f), Size: (%f, %f)\n", r.pos.x, r.pos.y, r.size.x, r.size.y);
}


static float rects_intersect(Rect r1, Rect r2) {
	if (r1.pos.x >= r2.pos.x + r2.size.x) return false; // r1 is to the right of r2
	if (r2.pos.x >= r1.pos.x + r1.size.x) return false; // r2 is to the right of r1
	if (r1.pos.y >= r2.pos.y + r2.size.y) return false; // r1 is above r2
	if (r2.pos.y >= r1.pos.y + r1.size.y) return false; // r2 is above r1
	return true;
}

// returns whether or not there is any of the clipped rectangle left
static bool rect_clip_to_rect(Rect *clipped, Rect clipper) {
	v2 start_pos = clipped->pos;
	clipped->pos.x = maxf(clipped->pos.x, clipper.pos.x);
	clipped->pos.y = maxf(clipped->pos.y, clipper.pos.y);
	clipped->size = v2_add(clipped->size, v2_sub(start_pos, clipped->pos));

	clipped->size.x = clampf(clipped->size.x, 0, clipper.pos.x + clipper.size.x - clipped->pos.x);
	clipped->size.y = clampf(clipped->size.y, 0, clipper.pos.y + clipper.size.y - clipped->pos.y);
	return clipped->size.x > 0 && clipped->size.y > 0;
}

// removes `amount` from all sides of r
static Rect rect_shrink(Rect r, float amount) {
	r.pos.x += amount;
	r.pos.y += amount;
	r.size.x -= 2 * amount;
	r.size.y -= 2 * amount;
	r.size.x = maxf(r.size.x, 0);
	r.size.y = maxf(r.size.y, 0);
	return r;
}

// adds `amount` to all sides of r
static Rect rect_grow(Rect r, float amount) {
	r.pos.x -= amount;
	r.pos.y -= amount;
	r.size.x += 2 * amount;
	r.size.y += 2 * amount;
	return r;
}

static v4 color_rgba_to_hsva(v4 rgba) {
	float R = rgba.x, G = rgba.y, B = rgba.z, A = rgba.w;
	float M = maxf(R, maxf(G, B));
	float m = minf(R, minf(G, B));
	float C = M - m;
	float H = 0;
	if (C == 0)
		H = 0;
	if (M == R)
		H = fmodf((G - B) / C, 6);
	else if (M == G)
		H = (B - R) / C + 2;
	else if (M == B)
		H = (R - G) / C + 4;
	H *= 60;
	float V = M;
	float S = V == 0 ? 0 : C / V;
	return V4(H, S, V, A);
}

static v4 color_hsva_to_rgba(v4 hsva) {
	float H = hsva.x, S = hsva.y, V = hsva.z, A = hsva.w;
	H /= 60;
	float C = S * V;
	float X = C * (1 - fabsf(fmodf(H, 2) - 1));
	float R, G, B;
	if (H <= 1)
		R=C, G=X, B=0;
	else if (H <= 2)
		R=X, G=C, B=0;
	else if (H <= 3)
		R=0, G=C, B=X;
	else if (H <= 4)
		R=0, G=X, B=C;
	else if (H <= 5)
		R=X, G=0, B=C;
	else
		R=C, G=0, B=X;
	
	float m = V-C;
	R += m;
	G += m;
	B += m;
	return V4(R, G, B, A);
}

static u32 color_interpolate(float x, u32 color1, u32 color2) {
	x = x * x * (3 - 2*x); // hermite interpolation

	v4 c1 = rgba_u32_to_v4(color1), c2 = rgba_u32_to_v4(color2);
	// to make it interpolate more nicely, convert to hsv, interpolate in that space, then convert back
	c1 = color_rgba_to_hsva(c1);
	c2 = color_rgba_to_hsva(c2);
	float h1 = c1.x, s1 = c1.y, v_1 = c1.z, a1 = c1.w;
	float h2 = c2.x, s2 = c2.y, v_2 = c2.z, a2 = c2.w; // avoid shadowing with v2
	
	float s_out = lerpf(x, s1, s2);
	float v_out = lerpf(x, v_1, v_2);
	float a_out = lerpf(x, a1, a2);
	
	float h_out;
	// because hue is on a circle, we need to make sure we take the shorter route around the circle
	if (fabsf(h1 - h2) < 180) {
		h_out = lerpf(x, h1, h2);
	} else if (h1 > h2) {
		h_out = lerpf(x, h1, h2 + 360);
	} else {
		h_out = lerpf(x, h1 + 360, h2);
	}
	h_out = fmodf(h_out, 360);
	
	v4 c_out = V4(h_out, s_out, v_out, a_out);
	c_out = color_hsva_to_rgba(c_out);
	return v4_to_rgba_u32(c_out);
}

// generate a random color with normally distributed red, green, and blue values
// the alpha value is taken directly from the mean
static u32 color_rand_gauss_seed(u64 *seed, u32 mean, float variation) {
	float channels[4];
	rgba_u32_to_floats(mean, channels);
	channels[0] += rand_gauss_seed(seed) * variation;
	channels[1] += rand_gauss_seed(seed) * variation;
	channels[2] += rand_gauss_seed(seed) * variation;
	channels[0] = clampf(channels[0], 0, 1);
	channels[1] = clampf(channels[1], 0, 1);
	channels[2] = clampf(channels[2], 0, 1);
	return floats_to_rgba_u32(channels);
}

static float barycentric_interpolation_2d(v2 p, v2 a, float a_v, v2 b, float b_v, v2 c, float c_v) {
	float w_a = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y))
		/ ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	float w_b = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y))
		/ ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	float w_c = 1 - w_a - w_b;
	assert(w_a >= 0 && w_a <= 1);
	assert(w_b >= 0 && w_b <= 1);
	assert(w_c >= 0 && w_c <= 1);
	return a_v * w_a + b_v * w_b + c_v * w_c;
}
