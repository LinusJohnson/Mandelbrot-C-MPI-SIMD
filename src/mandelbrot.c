#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include "mandelbrot.h"

void init_settings(int argc, char **argv, Settings *sett) {
	int rank, size;
	int SIMD = 1;
	int test = 1;
	int num_frames = 1;
	int frame_rate = 30;
	Frame_Res full_frame_res = {2048,2048};
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int cmd = 0;
	while ((cmd = getopt(argc, argv, "nt:v:f:h")) != -1) {
		switch (cmd) {
		case 'n':
			if(rank == ROOT) {printf("Running in non-SIMD mode.\n");}
			SIMD = 0;
			break;
		case 't':
			if(rank == ROOT) {printf("Running in performance test mode with %i readings \n", atoi(optarg));}
			test = atoi(optarg);
			if (test <= 0){
				if(rank == ROOT) {fprintf(stderr, "Option -%c requires a non-zero integer as argument.\n", optopt);}
				abort();
			}
			break;
		case 'v':
			if(rank == ROOT) {printf("Creating a video of %i frames \n", atoi(optarg));}
			num_frames = atoi(optarg);
			if (num_frames <= 0){
				if(rank == ROOT) {fprintf(stderr, "Option -%c requires a non-zero integer as argument.\n", optopt);}
				abort();
			}
			break;
		case 'f':
			if(rank == ROOT) {printf("Setting frame rate to %i \n", atoi(optarg));}
			frame_rate = atoi(optarg);
			if (frame_rate <= 0){
				if(rank == ROOT) {fprintf(stderr, "Option -%c requires a non-zero integer as argument.\n", optopt);}
				abort();
			}
			break;
		case 'r':
			if(rank == ROOT) {printf("Setting resolution to %i x %i \n", atoi(optarg), atoi(optarg));}
			full_frame_res.x = atoi(optarg);
			full_frame_res.y = atoi(optarg);
			if (frame_rate <= 0){
				if(rank == ROOT) {fprintf(stderr, "Option -%c requires a non-zero integer as argument.\n", optopt);}
				abort();
			}
			break;
		case 'h':
			if(rank == ROOT){ printf("\nThis program renders images and fractal zoom animation of the Mandelbrot set. \
										\nAuthor: Linus Johnson, ljohnson@kth.se \
										\nUse option: \n\t-n to run in non-SIMD mode (SIMD mode on by default) \
										\n\t-v n to render a video with n frames \
										\n\t-f n to specify the frame rate (default 30) \
										\n\t-r n to set the resolution to n x n (default 2048 x 2048) \
										\n\t-t n to run in performance test mode with n tests (must have been compiled with TIME_FRAMES set to 1 in definitions.h)\n\n");}
			MPI_Finalize();
			exit(0);
			break;
		case '?':
			if (optopt == 't') {
				if(rank == ROOT) {fprintf(stderr, "Option -%c requires an argument.\n", optopt);}
			} else if (isprint(optopt)) {
				if(rank == ROOT) {fprintf(stderr, "Unknown option `-%c'.\n", optopt);}
			} else {
				if(rank == ROOT) {fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);}
			}
	     	abort();
	     	break;
		default:
			break;
		}
	}
	if(SIMD & (rank == ROOT)) printf("Running in SIMD mode, with %s instructions and bit-width %i \n", SIMD_MODE, SIMD_WIDTH);

	const int sub_res_x = full_frame_res.x;
	if(full_frame_res.y % size != 0) { 
		fprintf(stderr, "The resolution must be divisible by the number of threads used\n");
		abort();
	}
	const int sub_res_y = full_frame_res.y / size;
	const Frame_Res sub_res = {sub_res_x, sub_res_y};
	const float b = 2;
	const float N = 255;
	const float dx = 2*b / full_frame_res.x;
	const float dy = 2*b / full_frame_res.y;
	const gen_pixel_func func = SIMD ? &gen_pixels_simd_intrinsics : &gen_pixels_no_simd;

	Colormap cmap;
	if (readColormap(&cmap, "viridis")) {
		fprintf(stderr, "Failed to read colormap\n");
		abort();
	}
	Frame_Settings f_s;
	memcpy(f_s.cmap, cmap, sizeof(f_s.cmap));
	f_s.zreal = 0; f_s.zimag = 0; f_s.zoom = 0;
	if ((num_frames < frame_rate) & (rank == ROOT))
		printf("To few frames setting frame rate to 1\n");
	frame_rate = num_frames < frame_rate ? 1 : frame_rate;
	const Video_Settings v_s = num_frames > 1 ? \
			(Video_Settings) {num_frames, frame_rate, num_frames / frame_rate, zoom_func_exp} : \
			(Video_Settings) {num_frames, 0, 0, zoom_func_exp};
	const Mandelbrot_Settings m_s = (Mandelbrot_Settings) {rank, size, full_frame_res, sub_res, b, N, dx, dy, func, test};
	sett->frame_settings = f_s;
	sett->video_settings = v_s;
	sett->mandelbrot_settings = m_s;
}

float fill_sub_frame(Frame_Data *f, const Settings *s) {
	TIME_START();
	gen_pixels(f, s);
	return TIME_END();
}

float fill_frame(Frame_Data *f_d, const Frame_Data *s_f_d, const Settings *sett) {
	// recieve subframes (with s.s_f_res.x offset since ROOT writes directly to the frame)
	// for ROOT s_f_d == NULL 
	int data_size = sett->mandelbrot_settings.s_f_res.x*sett->mandelbrot_settings.s_f_res.y;
	int frame_size = sett->mandelbrot_settings.f_res.x*sett->mandelbrot_settings.f_res.y;
	if (data_size == frame_size) return 0; // f_d already contains all data
	int sendcount = data_size;
	const int size = sett->mandelbrot_settings.size;
	const int rank = sett->mandelbrot_settings.rank;
	Frame_Data sendbuf = NULL;
	Frame_Data recvbuf = NULL;
	int *displs = NULL;;
	int *recvcounts = NULL;
	if (rank == ROOT){
		recvcounts = malloc(size*sizeof(int));
		displs = malloc(size*sizeof(int));
		for (int i = 0; i < size; ++i) {
			displs[i] = i*data_size;
			recvcounts[i] = data_size;
		}
		recvcounts[0] = 0; // root does not send any data
		recvbuf = &((*f_d)[0]);
		sendcount=0;
	} else {
		sendbuf = &((*s_f_d)[0]);
	}
	TIME_START();
	MPI_Gatherv(sendbuf, sendcount, MPI_UNSIGNED, 
				recvbuf, recvcounts, displs, MPI_UNSIGNED,  
				ROOT, MPI_COMM_WORLD); 
	float t = TIME_END();
	if(rank==ROOT) {
		free(recvcounts);
		free(displs);
	}
	return t;
}

void get_frame(Frame *f, const Settings *sett) {
	const int rank = sett->mandelbrot_settings.rank;
	Frame_Data *f_d = NULL; // NULL if not ROOT
	Frame_Data *sub_f_d = NULL; // NULL if ROOT
	if (rank == ROOT) {
		f_d = malloc(sizeof(Frame_Data));
		malloc2d(f_d, sett->mandelbrot_settings.f_res.x, 
					  sett->mandelbrot_settings.f_res.y);
	} else {
		sub_f_d = malloc(sizeof(Frame_Data));
		malloc2d(sub_f_d, sett->mandelbrot_settings.s_f_res.x, 
						  sett->mandelbrot_settings.s_f_res.y);
	}
	#if TIME_FRAMES
	const float calc_time = fill_sub_frame((rank == ROOT) ?  f_d : sub_f_d, sett);
	const float comm_time = fill_frame(f_d, sub_f_d, sett);
	Frame_Time f_t = {calc_time, comm_time};
	if (rank == ROOT) {
		f->data = *f_d;
		f->t = f_t;
	}
	#else
	fill_sub_frame((rank == ROOT) ?  f_d : sub_f_d, sett);
	fill_frame(f_d, sub_f_d, sett);
	if (rank == ROOT) f->data = *f_d;
	#endif
	if (rank != ROOT) {
		free_frame_data(sub_f_d);
	}
}

void finalize(Frame *frame_pointer, Settings *sett) {
	MPI_Finalize();
	if (sett->mandelbrot_settings.rank == ROOT){
		if (sett->video_settings.num_frames > 1)
			free(frame_pointer); // frame data already freed
		else
			free_frame(frame_pointer);
	}
}

void free_frame(Frame *frame_pointer) {
	free_frame_data(&(frame_pointer->data));
	free(frame_pointer);
}

void gen_pixels(Frame_Data *f_d, const Settings *sett) {
	sett->mandelbrot_settings.func(*f_d, sett);
}

float zoom_func_exp(int frame_num) {
	return 4. + exp((frame_num / 60.));
}