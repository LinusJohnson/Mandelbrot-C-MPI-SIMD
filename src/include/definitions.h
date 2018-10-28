#ifndef DEFINITIONS_
#define DEFINITIONS_
#include <stdint.h>
#define TIME_FRAMES 0
#define ROOT 0
typedef float c_t;
typedef uint32_t Datatype;
typedef Datatype *Frame_Data;
typedef float Colormap[256][3];

typedef struct Frame_Res {
	int x;
	int y;
} Frame_Res;

typedef struct Frame_Settings {
	float zreal;
	float zimag;
	float zoom;
	Colormap cmap;
} Frame_Settings;

typedef float(*zoom_scale_func)(int);

typedef struct Video_Settings {
	int num_frames;
	int frame_rate;
	float duration;
	zoom_scale_func zoom_func;
} Video_Settings;

struct Settings;

typedef struct Mandelbrot_Settings {
	int rank;
	int size;
	Frame_Res f_res; // full frame resolution
	Frame_Res s_f_res; // sub frame resolution
	float b;
	int N;
	float dx;
	float dy;
	void (*func)(Frame_Data, const struct Settings*);
	int test;
} Mandelbrot_Settings;

typedef struct Settings {
	Frame_Settings frame_settings;
	Video_Settings video_settings;
	Mandelbrot_Settings mandelbrot_settings;
} Settings;

typedef void(*gen_pixel_func)(Frame_Data, const Settings*);



#define LOWEST_FLOAT_PREC (4*__CHAR_BIT__)
#if defined(__AVX512F__)
#define SIMD_WIDTH 512
#define SIMD_MODE "AVX512F"
#elif defined(__AVX2__)
#define SIMD_WIDTH 256
#define SIMD_MODE "AVX2"
#elif defined(__AVX__)
#define SIMD_WIDTH 256
#define SIMD_MODE "AVX"
#elif defined(__SSE3__)
#define SIMD_WIDTH 128
#define SIMD_MODE "SSE3"
#elif defined(__SSE2__)
#define SIMD_WIDTH 128
#define SIMD_MODE "SSE2"
#elif defined(__SSE__)
#define SIMD_WIDTH 128
#define SIMD_MODE "SSE"
#else
#define SIMD_WIDTH 32
#define SIMD_MODE "Unknown"
#endif
#define SIMD_N (SIMD_WIDTH / LOWEST_FLOAT_PREC)

#if TIME_FRAMES
typedef struct Frame_Time {
	float calc_t;
	float comm_t;
} Frame_Time;
#endif

#if TIME_FRAMES
typedef struct Frame {
	Frame_Data data;
	Frame_Time t;
} Frame;
#else
typedef struct Frame {
	Frame_Data data;
} Frame;
#endif

#endif