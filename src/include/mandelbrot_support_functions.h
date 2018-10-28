#ifndef MANDELBROT_SUPPORT_FUNCTIONS_
#define MANDELBROT_SUPPORT_FUNCTIONS_
#include <stdio.h>
#include <png.h>
#include "definitions.h"
#include "mandelbrot.h"
#include "mandelbrot_functions.h"

int writeImage(char* filename, const Settings *settings, Frame_Data buffer, char* title);
void malloc2d(Frame_Data *array, int n, int m);
void free_frame_data(Frame_Data *array);
void print_to_file(Frame_Data f_d, const Settings *settings);
int render_video(const char *filename, Frame *frames, Settings *settings);
int readColormap(Colormap *cmap, char *filename);

#endif
