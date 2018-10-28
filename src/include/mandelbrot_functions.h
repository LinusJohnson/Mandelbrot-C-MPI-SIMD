#ifndef MANDELBROT_FUNCIONS_
#define MANDELBROT_FUNCIONS_
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <immintrin.h>
#include "definitions.h"
#include "mandelbrot.h"

void gen_pixels_no_simd(Frame_Data __restrict f_d, const Settings *sett);
void gen_pixels_simd(Frame_Data __restrict f_d, const Settings *sett);
void gen_pixels_simd_intrinsics(Frame_Data __restrict f_d, const Settings *sett);

#if SIMD_WIDTH==256
#define _mul_ps(a,b) _mm256_mul_ps(a,b)
#define _add_ps(a,b) _mm256_add_ps(a,b)
#define _sub_ps(a,b) _mm256_sub_ps(a,b)
#define _or_ps(a,b) _mm256_or_ps(a,b)
#define _set1_ps(a) _mm256_set1_ps(a)
#define _cmplt_ps(a,b) _mm256_cmp_ps(a,b,_CMP_LT_OS)
#define _cmpeq_ps(a,b) _mm256_cmp_ps(a,b,_CMP_EQ_OS)
#define _cmpgt_ps(a,b) _mm256_cmp_ps(a,b,_CMP_GT_OS)
#define _or_si(a,b) _mm256_or_si256(a,b)
#define _and_si(a,b) _mm256_and_si256(a,b)
#define _and_ps(a,b) _mm256_and_ps(a,b)
#define _andnot_si(a,b) _mm256_andnot_si256(a,b)
#define _andnot_ps(a,b) _mm256_andnot_ps(a,b)
#define _extract_ps(a,b) _mm256_extract_ps(a,b)
#define _store_si(a,b) _mm256_store_si256(a,b)
#define _set_ps(a,b,c,d,e,f,g,h) _mm256_set_ps(a,b,c,d,e,f,g,h)
#define _castps_si(a) _mm256_castps_si256(a)
#define _xor_ps(a,b) _mm256_xor_ps(a,b)

#define _setzero_ps() _mm256_setzero_ps()
#define _cvtss_f32(a) _mm256_cvtss_f32(a)
#if defined(__AVX2__)
#define _cvtsi_si32(a) _mm256_cvtsi256_si32(a)
#define _movemask_epi8(a) _mm256_movemask_epi8(a)
#else
#define _cvtsi_si32(a) (int) _mm256_cvtss_f32(_mm256_castsi256_ps(a))
#define _movemask_epi8(a) _mm256_movemask_ps(_mm256_castsi256_ps(a))
#endif
#define _shuffle_ps(a,b,c) _mm256_shuffle_ps(a,b,c)
#define _set_ps_x(x) _set_ps(x+7,x+6,x+5,x+4,x+3,x+2,x+1,x)

#define _storeu_si(a,b) _mm256_storeu_si256(a,b)
#define _storeu_ps(a,b) _mm256_storeu_ps(a,b)

//#if LOWEST_FLOAT_PREC==32
#define _shuffle_epi32(a,b) _mm256_shuffle_epi32(a,b)
#define _add_epi32(a,b) _mm256_add_epi32(a,b)
#define _cvtps_epi32(a) _mm256_cvtps_epi32(a)
#define _set1_epi32(a) _mm256_set1_epi32(a)
#define _cmpgt_epi32(a,b) _mm256_cmpgt_epi32(a,b)
#define _extract_epi32(a,b) _mm256_extract_epi32(a,b)

typedef __m256i __mi;
typedef __m256 __m;

static inline void __print_epi32(char * str, __mi in) {
    uint32_t *val = (uint32_t*) &in;
    printf("%s %i %i %i %i %i %i %i %i \n", str, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
}

static inline void __print_ps(char * str, __m in) {
    float *val = (float*) &in;
    printf("%s %f %f %f %f %f %f %f %f \n", str, val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7]);
}

#else

#define _mul_ps(a,b) _mm_mul_ps(a,b)
#define _add_ps(a,b) _mm_add_ps(a,b)
#define _add_epi32(a,b) _mm_add_epi32(a,b)
#define _sub_ps(a,b) _mm_sub_ps(a,b)
#define _cvtps_epi32(a) _mm_cvtps_epi32(a)
#define _or_ps(a,b) _mm_or_ps(a,b)
#define _set1_ps(a) _mm_set1_ps(a)
#define _cmplt_ps(a,b) _mm_cmplt_ps(a,b)
#define _cmpeq_ps(a,b) _mm_cmpeq_ps(a,b)
#define _cmpgt_ps(a,b) _mm_cmpgt_ps(a,b)
#define _or_si(a,b) _mm_or_si128(a,b)
#define _and_si(a,b) _mm_and_si128(a,b)
#define _and_ps(a,b) _mm_and_ps(a,b)
#define _set1_epi32(a) _mm_set1_epi32(a)
#define _andnot_si(a,b) _mm_andnot_si128(a,b)
#define _andnot_ps(a,b) _mm_andnot_ps(a,b)
#define _extract_epi32(a,b) _mm_extract_epi32(a,b)
#define _extract_ps(a,b) _mm_extract_ps(a,b)
#define _cmpgt_epi32(a,b) _mm_cmpgt_epi32(a,b)
#define _store_si(a,b) _mm_store_si128(a,b)
#define _set_ps(a,b,c,d) _mm_set_ps(a,b,c,d)
#define _castps_si(a) _mm_castps_si128(a)
#define _xor_ps(a,b) _mm_xor_ps(a,b)
#define _movemask_epi8(a) _mm_movemask_epi8(a)
#define _setzero_ps() _mm_setzero_ps()
#define _cvtss_f32() _mm_cvtss_f32()
#define _cvtsi_si32(a) _mm_cvtsi128_si32(a)
#define _shuffle_ps(a,b,c) _mm_shuffle_ps(a,b,c)
#define _shuffle_epi32(a,b) _mm_shuffle_epi32(a,b)
#define _set_ps_x(x) _set_ps(x+3,x+2,x+1,x)
#define _storeu_si(a,b) _mm_storeu_si128(a,b)
#define _storeu_ps(a,b) _mm_storeu_ps(a,b)

typedef __m128i __mi;
typedef __m128 __m;

typedef void (*__print_handler)(char*, __m);

void __print_epi32(char * str, __mi in) {
    uint32_t *val = (uint32_t*) &in;
    printf("%s %i %i %i %i \n", str, val[0], val[1], val[2], val[3]);
}

void __print_ps(char * str, __m in) {
    float *val = (float*) &in;
    printf("%s %f %f %f %f \n", str, val[0], val[1], val[2], val[3]);
}

#endif

#define __print(a,b) _Generic((b), __m: __print_ps, __mi: __print_epi32)(a,b)

#define _not_ps(x) _xor_ps(x, _cmpeq_ps(_setzero_ps(), _setzero_ps()))

// (zreal_sq[n] - 0.5*zreal[n] + 0.0625 + zimag_sq)
#define Q() _add_ps(_add_ps(_mul_ps(_set1_ps(-0.5), zreal), _add_ps(zreal_sq, _set1_ps(0.0625))), zimag_sq) 


//#pragma omp declare simd uniform(dreal, b,  N, SIMD_N) notinbranch simdlen(4)
static inline void cal_pixel_simd(Frame_Data __restrict count, const c_t dreal, const c_t *__restrict dimag, 
						const int b, const int N) {
	c_t zreal[SIMD_N];
	c_t zimag[SIMD_N];
	#pragma omp simd
	for(int n=0; n < SIMD_N; n++){zreal[n] = dreal;}
	memcpy(zimag,dimag,SIMD_N*sizeof(c_t));
	c_t zreal_prev[SIMD_N];
	c_t zimag_prev[SIMD_N];
	c_t z_abs[SIMD_N];
	c_t zreal_sq[SIMD_N];
	c_t zimag_sq[SIMD_N];
	c_t zri[SIMD_N];
	int run[SIMD_N];
	// bulb checking (see https://en.wikipedia.org/wiki/Mandelbrot_set#Cardioid_/_bulb_checking)
	// ~4.5 times speedup
	#pragma omp simd
	for(int n=0; n < SIMD_N; n++){
		zreal_sq[n] = zreal[n]*zreal[n];
		zimag_sq[n] = zimag[n]*zimag[n];
		run[n] = !(((zreal_sq[n] - 0.5*zreal[n] + 0.0625) + zimag_sq[n])*(((zreal_sq[n] - 0.5*zreal[n] + 0.0625) + zimag_sq[n]) + (zreal[n]-0.25)) < 0.25*zimag_sq[n]) || 
				 (zreal_sq[n]+2*zreal[n]+1+zimag_sq[n] < 0.0625);
		count[n] = run[n] ? 0 : N;
	}
	int done = 0;
	#pragma omp simd reduction(|:done)
	for(int n=0; n<SIMD_N; n++){done|=run[n];}
	done=!done;
	while(!done){
		#pragma omp simd
		for(int n=0; n<SIMD_N; n++) {  
			zreal_prev[n] = zreal[n];
			zimag_prev[n] = zimag[n];
			zreal_sq[n] = zreal[n]*zreal[n];
			zimag_sq[n] = zimag[n]*zimag[n];
			z_abs[n] = zreal_sq[n] + zimag_sq[n];
			zri[n] = zreal[n]*zimag[n];
			zreal[n] = zreal_sq[n] - zimag_sq[n] + dreal;
			zimag[n] = 2*zri[n] + dimag[n];
			count[n] += run[n];
			// periodicity checking (see https://en.wikipedia.org/wiki/Mandelbrot_set#Periodicity_checking)
			// ~1.07 times speedup
			count[n] = (zreal[n] == zreal_prev[n] && zimag[n] == zimag_prev[n]) ? N : count[n];
			run[n] &= (count[n] < N) & (z_abs[n] < b);
		}
		#pragma omp simd reduction(|:done)
		for(int n=0; n<SIMD_N; n++){done|=run[n];}
		done=!done;
	}
}

static inline void cal_pixel_simd_intrinsics(__mi *__restrict f_d, const __m dreal, const __m dimag, 
											 const int b, const int N) {
	__m count = _set1_ps(0);
	__m zreal = dreal;
	__m zimag = dimag;
	__m zreal_prev;
	__m zimag_prev;
	__m z_abs;
	__m zreal_sq;
	__m zimag_sq;
	__m zri;
	__m run;
	// bulb checking (see https://en.wikipedia.org/wiki/Mandelbrot_set#Cardioid_/_bulb_checking)
	// ~4.5 times speedup
	//zreal_sq[n] = zreal[n]*zreal[n];
	zreal_sq = _mul_ps(zreal,zreal);
	//zimag_sq[n] = zimag[n]*zimag[n];
	zimag_sq = _mul_ps(zimag,zimag);
	//run[n] = !(((zreal_sq[n] - 0.5*zreal[n] + 0.0625) + zimag_sq[n])*
	// 			(((zreal_sq[n] - 0.5*zreal[n] + 0.0625) + zimag_sq[n]) + (zreal[n]-0.25)) < 0.25*zimag_sq[n]) || 
	//			 (zreal_sq[n]+2*zreal[n]+1+zimag_sq[n] < 0.0625);
	run = _not_ps(_or_ps(_cmplt_ps(_mul_ps(Q(), _add_ps(Q(), _add_ps(zreal, _set1_ps(-0.25)))), 
								   _mul_ps(_set1_ps(0.25), zimag_sq)),
					 _cmplt_ps(_add_ps(_add_ps(_mul_ps(_set1_ps(2), zreal), zreal_sq), _add_ps(_set1_ps(1), zimag_sq)), 
								   _set1_ps(0.0625))));
	//count[n] = run[n] ? 0 : N;
	count = _or_ps(_and_ps(count, run), _andnot_ps(run, _set1_ps(N)));
	int done = 0;
	u_int16_t test = _movemask_epi8(_castps_si(run));
	done = !((test != 0) | (test == 0xffff));
	while(!done){
		// zreal_prev[n] = zreal[n];
		zreal_prev = zreal;
		// zimag_prev[n] = zimag[n];
		zimag_prev = zimag;
		// zreal_sq[n] = zreal[n]*zreal[n];
		zreal_sq = _mul_ps(zreal,zreal);
		// zimag_sq[n] = zimag[n]*zimag[n];
		zimag_sq = _mul_ps(zimag,zimag);
		// z_abs[n] = zreal_sq[n] + zimag_sq[n];
		z_abs = _add_ps(zreal_sq, zimag_sq);
		// zri[n] = zreal[n]*zimag[n];
		zri = _mul_ps(zreal,zimag);
		// zreal[n] = zreal_sq[n] - zimag_sq[n] + dreal;
		zreal = _add_ps(_sub_ps(zreal_sq, zimag_sq), dreal);
		// zimag[n] = 2*zri[n] + dimag[n];
		zimag = _add_ps(_mul_ps(_set1_ps(2), zri), dimag);
		// count[n] += run[n];
		count = _add_ps(count, _and_ps(run, _set1_ps(1)));
		// count[n] = (zreal[n] == zreal_prev[n] && zimag[n] == zimag_prev[n]) ? N : count[n];
		__m period_mask = _and_ps(_cmpeq_ps(zreal, zreal_prev), _cmpeq_ps(zimag, zimag_prev));
		count = _or_ps(_and_ps(_set1_ps(N), period_mask), _andnot_ps(period_mask, count));
		// run[n] &= (count[n] < N) & (z_abs[n] < b);
		__m z_abs_mask = _cmplt_ps(z_abs, _set1_ps(b));
		__m count_mask = _cmplt_ps(count, _set1_ps(N));
		run = _and_ps(_and_ps(z_abs_mask, count_mask), run);
		test =  _movemask_epi8(_castps_si(run));
		done = !((test != 0) | (test == 0xffff));
	}
	_storeu_si(f_d, _cvtps_epi32(count));
}

static inline Datatype cal_pixel_no_simd(c_t dreal, c_t dimag, int b, int N) {
	Datatype count = 0;
	c_t zreal = dreal;
	c_t zimag = dimag;
	c_t zreal_prev;
	c_t zimag_prev;
	c_t z_abs;
	int run = 1;
	c_t q = (zreal*zreal - 0.5*zreal + 0.0625) + zimag*zimag;
	// bulb checking (see https://en.wikipedia.org/wiki/Mandelbrot_set#Cardioid_/_bulb_checking)
	// ~4.5 times speedup
	if ((q*(q + (zreal-0.25)) < 0.25*zimag*zimag) || (zreal*zreal+2*zreal+1+zimag*zimag < 0.0625)) {
		run = 0; count = N;
	}
	while(run){
		zreal_prev = zreal;
		zimag_prev = zimag;
		c_t zreal_sq = zreal*zreal;
		c_t zimag_sq = zimag*zimag;
		z_abs = zreal_sq + zimag_sq;
		c_t zri = zreal*zimag;
		zreal = zreal_sq - zimag_sq + dreal;
		zimag = 2*zri + dimag;
		count++;
		// periodicity checking (see https://en.wikipedia.org/wiki/Mandelbrot_set#Periodicity_checking)
		// ~1.07 times speedup
		count = (zreal == zreal_prev && zimag == zimag_prev) ? N : count;
		run = (count < N) & (z_abs < b);
	}
	return count;
}

#endif