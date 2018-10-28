#include <stdlib.h>
#include <stdio.h>

#include "mandelbrot_main.h"
#include "mandelbrot_support_functions.h"
#include "mandelbrot.h"

int main(int argc, char **argv) {
	Settings s;
	init_settings(argc, argv, &s);
	const int rank = s.mandelbrot_settings.rank;
	const int test = s.mandelbrot_settings.test;
	const int num_frames = s.video_settings.num_frames;
	Frame *frame = NULL;
	if (rank == ROOT) frame = malloc(sizeof(frame));
	#if TIME_FRAMES
	float calc_t = 0;
	float comm_t = 0;
	#endif
	s.frame_settings.zreal = -0.743639266077433;
	s.frame_settings.zimag = 0.131824786875559;
	s.frame_settings.zoom = 1;
	if ((num_frames > 1) & (test > 1)) {
		fprintf(stderr, "Performance testing mode doesn't support video rendering");
		abort();
	}
	#if !TIME_FRAMES
	if (test > 1) {
		fprintf(stderr, "Can't run performance testing, program hasn't been compiled for timing!");
		abort();
	}
	#endif
	if (num_frames > 1) {
		if (rank == ROOT) {
			printf("Rendering video...\n");
			render_video("mandelbrot.mp4", frame, &s);
			printf("Rendering completed!\n");
			printf("Animation saved as mandelbrot.mp4\n");
		} else {
			for (int i=0; i < num_frames; i++){
				// the other processes continue to get frames 
				// and send them to ROOT through Gatherv while 
				// ROOT goes through the video rendering  
				get_frame(frame, &s);
				s.frame_settings.zoom = s.video_settings.zoom_func(i);
			} 
		}
	} else {
		for (int i=0; i<test; i++){
			if (rank == ROOT) printf("Rendering image...\n");
			get_frame(frame, &s);
			if (rank == ROOT) printf("Rendering completed!\n");
			#if TIME_FRAMES
			if (rank == ROOT) {
				calc_t+=frame->t.calc_t;
				comm_t+=frame->t.comm_t;
			}
			#endif
		}
	}
	#if TIME_FRAMES
	calc_t/=test;
	comm_t/=test;
	#endif
	if ((rank == ROOT) & (num_frames == 1)) {
		print_to_file(frame->data, &s);
		printf("Data saved to text file mandelbrot.txt\n");
		writeImage("mandelbrot.png", &s, frame->data, "mandelbrot");
		printf("Image saved as mandelbrot.png\n");
		#if TIME_FRAMES
		if(test>1) {
			printf("Average calculation time: %f s\n", calc_t);
			printf("Average communication time: %f s\n", comm_t);
			printf("Average total Execution time: %f s\n", calc_t + comm_t);
		} else { 
			printf("Calculation time: %f s\n", calc_t);
			printf("Communication time: %f s\n", comm_t);
			printf("Total Execution time: %f s\n", calc_t + comm_t);
		}
		#endif
	}
	finalize(frame,&s);
	return 0;
}