#include <math.h>
#include <limits.h>
#include "mandelbrot.h"
#include <mpi.h>

void gen_pixels_no_simd(Frame_Data __restrict f_d, const Settings *sett) {
	const Mandelbrot_Settings *s = &(sett->mandelbrot_settings);
	const Frame_Settings f_sett = sett->frame_settings;
	int xoff = 0;
	int yoff = s->rank*s->s_f_res.y;
	int x,y;
	#pragma omp parallel for ordered schedule(dynamic)
	for(y=0; y < s->s_f_res.y; y++) {
		float dimag = (y+yoff-s->f_res.y/2.0)*s->dy/f_sett.zoom + f_sett.zimag;
		for(x=0; x < s->s_f_res.x; x++) {
			float dreal = (x+xoff-s->f_res.x/2.0)*s->dx/f_sett.zoom + f_sett.zreal;
			f_d[x + y*s->f_res.x] =  cal_pixel_no_simd(dreal, dimag, s->b, s->N);
		}
	}
}

void gen_pixels_simd(Frame_Data __restrict f_d, const Settings *sett) {
	const Mandelbrot_Settings *s = &(sett->mandelbrot_settings);
	const Frame_Settings f_sett = sett->frame_settings;
	int xoff = 0 ;
	int yoff = s->rank*s->s_f_res.y;
	int x,y;

	c_t dimag[s->s_f_res.y];
	#pragma omp simd
	for (y=0; y < s->s_f_res.y; y++) {dimag[y] = (y+yoff-s->f_res.y/2)*s->dy/f_sett.zoom - f_sett.zimag;}

	c_t dreal[s->s_f_res.x];
	#pragma omp simd
	for (x=0; x < s->s_f_res.x; x++) {dreal[x] = (x+xoff-s->f_res.x/2)*s->dx/f_sett.zoom - f_sett.zreal;}

	#pragma omp parallel for ordered collapse(2) schedule(dynamic)
	for(y=0; y < s->s_f_res.y; y++) {
		for(x=0; x < s->s_f_res.x; x+=SIMD_N) {
			Datatype pix[SIMD_N];
			cal_pixel_simd(pix, dreal[x], &dimag[y], s->b, s->N);
			#pragma omp simd
			for (int n=0; n < SIMD_N; n++) {
				f_d[x + y*s->f_res.x + n] = pix[n];
			}
		}
	}
}

void gen_pixels_simd_intrinsics(Frame_Data __restrict f_d, const Settings *sett) {
	const Mandelbrot_Settings *s = &(sett->mandelbrot_settings);
	const Frame_Settings f_sett = sett->frame_settings;
	int xoff = 0;
	int yoff = s->rank*s->s_f_res.y;
	int x,y;
	//#pragma omp parallel for ordered schedule(dynamic)
	for(y=0; y < s->s_f_res.y; y++) {
		__m dimag = _set1_ps((y + yoff - s->f_res.y / 2.0)*s->dy / f_sett.zoom + f_sett.zimag);
		for(x=0; x < s->s_f_res.x; x+=SIMD_N) {
			__m dreal = _add_ps(_mul_ps(_add_ps(_set_ps_x(x), _set1_ps(xoff - s->f_res.x / 2.0)), _set1_ps(s->dx / f_sett.zoom)), _set1_ps(f_sett.zreal));
			// will calculate the x edge values twice if X_RES is not divisible by SIMD_N
			// f_d is padded to handle the extra bytes at the last row
			cal_pixel_simd_intrinsics((__mi *)&f_d[x + y*s->f_res.x], dreal, dimag, s->b, s->N);
		}
	}
}
