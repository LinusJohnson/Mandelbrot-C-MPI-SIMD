#ifndef MANDELBROT_
#define MANDELBROT_
#include <mpi.h>
#include <math.h>
#include "definitions.h"
#include "mandelbrot_support_functions.h"

void init_settings(int argc, char **argv, Settings *sett);

#if TIME_FRAMES
#define TIME_START() MPI_Barrier(MPI_COMM_WORLD); float start = MPI_Wtime();
#else
#define TIME_START()
#endif 

#if TIME_FRAMES
#define TIME_END_INNER() MPI_Barrier(MPI_COMM_WORLD); float end = MPI_Wtime(); duration = end - start;
#else
#define TIME_END_INNER()
#endif

#define TIME_END() ({float duration = 0; TIME_END_INNER(); duration;})

void free_frame(Frame *frame_pointer);
void finalize(Frame *frame_pointer, Settings *sett);

float fill_sub_frame(Frame_Data *f, const Settings *sett);
float fill_frame(Frame_Data *f_d, const Frame_Data *s_f_d, const Settings *sett);
void get_frame(Frame *f, const Settings *sett);
void gen_pixels(Frame_Data *f_d, const Settings *sett);

float zoom_func_exp(int frame_num);

#endif