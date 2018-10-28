#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libavcodec/avcodec.h>
#include <libavutil/opt.h>
#include <libavutil/imgutils.h>
#include <libswscale/swscale.h>
#include <libavutil/avassert.h>
#include <libavutil/channel_layout.h>
#include <libavutil/opt.h>
#include <libavutil/mathematics.h>
#include <libavutil/timestamp.h>
#include <libavformat/avformat.h>

#include "mandelbrot_support_functions.h"


#define CMAP_PATH "./src/data/colormaps/%s.%s"
#define SCALE_FLAGS SWS_LANCZOS

typedef struct OutputStream {
    AVStream *st;
    AVCodecContext *enc;
    int64_t next_pts;
    int samples_count;
    AVFrame *frame;
    AVFrame *tmp_frame;
    float t, tincr, tincr2;
    struct SwsContext *sws_ctx;
} OutputStream;

int readCSVColormap(Colormap* cmap, char* filepath, char* filepath_data) {
   FILE *fp;
   if ((fp = fopen(filepath, "r"))==NULL) {
      return 1;
   }
   size_t count = 0;
   for (; count < 256; ++count) {
      if (fscanf(fp, "%f,%f,%f", &((*cmap)[count][0]), &((*cmap)[count][1]), &((*cmap)[count][2])) != 3) {
         if(feof(fp)){
            fprintf(stderr, "EOF reached prematurely");
         } else {
            fprintf(stderr, "Could not read file %s", filepath);
         }
         return 1;
      }
   }
   FILE *fp_data;
   if((fp_data=fopen(filepath_data, "wb"))==NULL) {
      fprintf(stderr, "Could not open file %s for writing\n", filepath_data);
   }

   if(fwrite(cmap, sizeof(Colormap), 1, fp_data) != 1){
      fprintf(stderr, "Could not write to file %s", filepath_data);
   }

   fclose(fp_data);
   fclose(fp);
   return 0;
}

int readDATAColormap(Colormap *cmap, char *filepath) {
   FILE *fp;
   if ((fp = fopen(filepath, "rb")) == NULL) {
      return 1;
   } else {
      if(fread(cmap, sizeof(Colormap), 1, fp) != 1) {
         if(feof(fp)){
            fprintf(stderr, "EOF reached prematurely");
         } else {
            fprintf(stderr, "Could not read file %s", filepath);
         }
         return 1;
      }
   }
   fclose(fp);
   return 0;
}

int readColormap(Colormap *cmap, char *filename) {
   char filepath_data[80];
   sprintf(filepath_data, CMAP_PATH, filename, "data" );
	if(readDATAColormap(cmap, filepath_data)){
      char filepath_csv[80];
      sprintf(filepath_csv, CMAP_PATH, filename, "csv" );
      if(readCSVColormap(cmap, filepath_csv, filepath_data)){
         fprintf(stderr, "Failed to read file %s and %s \n", filepath_data, filepath_csv);
         return 1;
      }
   }
   return 0;
}

void setRGB(png_byte *ptr, Datatype val, const Colormap *cmap) {
	ptr[0] = 255*(*cmap)[val][0];
	ptr[1] = 255*(*cmap)[val][1];
	ptr[2] = 255*(*cmap)[val][2];
}

void free_png_ptrs(FILE **fp, png_structp *png_ptr, png_infop *info_ptr, png_bytep *row){
   if (*fp != NULL) fclose(*fp);
   if (*info_ptr != NULL) png_free_data(*png_ptr, *info_ptr, PNG_FREE_ALL, -1);
   if (*png_ptr != NULL) png_destroy_write_struct(png_ptr, (png_infopp)NULL);
   if (*row != NULL) free(*row);
}

int writeImage(char* filename, const Settings *settings, Frame_Data buffer, char* title){
   const int width = settings->mandelbrot_settings.f_res.x;
   const int height = settings->mandelbrot_settings.f_res.y; 
   FILE *fp = NULL;
   png_structp png_ptr = NULL;
   png_infop info_ptr = NULL;
   png_bytep row = NULL;
	// Open file for writing (binary mode)
   if ((fp = fopen(filename, "wb"))==NULL) {
      fprintf(stderr, "Could not open file %s for writing\n", filename);
      free_png_ptrs(&fp, &png_ptr,&info_ptr,&row);
      return 1;
   }
        // Initialize write structure
   if ((png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))==NULL) {
      fprintf(stderr, "Could not allocate write struct\n");
      free_png_ptrs(&fp, &png_ptr,&info_ptr,&row);
      return 1;
   }

   // Initialize info structure
   if ((info_ptr = png_create_info_struct(png_ptr))==NULL) {
      fprintf(stderr, "Could not allocate info struct\n");
      free_png_ptrs(&fp, &png_ptr,&info_ptr,&row);
      return 1;
   }
      // Setup Exception handling
   if (setjmp(png_jmpbuf(png_ptr))) {
      fprintf(stderr, "Error during png creation\n");
      free_png_ptrs(&fp, &png_ptr,&info_ptr,&row);
      return 1;
   }

   png_init_io(png_ptr, fp);

   // Write header (8 bit colour depth)
   png_set_IHDR(png_ptr, info_ptr, width, height,
         8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
         PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   // Set title
   if (title != NULL) {
      png_text title_text;
      title_text.compression = PNG_TEXT_COMPRESSION_NONE;
      title_text.key = "Title";
      title_text.text = title;
      png_set_text(png_ptr, info_ptr, &title_text, 1);
   }

   png_write_info(png_ptr, info_ptr);

   // Allocate memory for one row (3 bytes per pixel - RGB)
   row = (png_bytep) malloc(3 * width * sizeof(png_byte));


   // Write image data
   int x, y;
   for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
         setRGB(&(row[x*3]), buffer[x + y*height], &(settings->frame_settings.cmap));
      }
      png_write_row(png_ptr, row);
   }

   // End write
   png_write_end(png_ptr, NULL); 
   free_png_ptrs(&fp, &png_ptr, &info_ptr, &row);  
   return 0;
}

void print_to_file(Frame_Data f_d, const Settings *settings) {
   const int width = settings->mandelbrot_settings.f_res.x;
   const int height = settings->mandelbrot_settings.f_res.y; 
   FILE *fp;   
   fp = fopen("mandelbrot.txt", "w");
   for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
         fprintf(fp, "%u ", f_d[x + y*height]);
      }
      fprintf(fp, "\n");
   }
   fclose(fp);
}

void malloc2d(Frame_Data *array, int n, int m) {
   int size = n * m;
   int x_padding = (SIMD_N - (n % SIMD_N)) % SIMD_N;
   *array = aligned_alloc(SIMD_N, (size + x_padding) * sizeof(Datatype));
}

void free_frame_data(Frame_Data *frame_data) {
    free(*frame_data);
}


// Mostly taken from the FFMPEG examples https://github.com/FFmpeg/FFmpeg/tree/master/doc/examples

static int write_video_frame(AVFormatContext *oc, OutputStream *ost, Frame *frame_pointer, Settings *settings) {
   int ret;
   AVCodecContext *c;
   AVFrame *av_frame;
   int got_packet = 0;
   AVPacket pkt = { 0 };
   c = ost->enc;
   
   if (ost->next_pts >= settings->video_settings.num_frames)
      return 1;
   /* check if we want to generate more frames */
   if (av_compare_ts(ost->next_pts, c->time_base, settings->video_settings.duration, (AVRational){ 1, 1 }) >= 0)
      return 1;

   /* when we pass a frame to the encoder, it may keep a reference to it
   * internally; make sure we do not overwrite it here */
   if (av_frame_make_writable(ost->frame) < 0) exit(1);

   ost->sws_ctx = sws_getContext(c->width, c->height,
                                 AV_PIX_FMT_RGB24, c->width, c->height,
                                 AV_PIX_FMT_YUV420P, SCALE_FLAGS, NULL, NULL, NULL);
   uint8_t *RGB_data = NULL;
   get_frame(frame_pointer, settings);

   RGB_data = malloc(3 * c->width * c->height * sizeof(uint8_t));
   int x,y;
   for (y = 0; y < c->height; y++) {
      for (x = 0; x < c->width; x++) {
         setRGB(&(RGB_data[3*(y*c->width + x)]), frame_pointer->data[x + y*c->height], &(settings->frame_settings.cmap));
      }
   }

   int inLinesize[1] = { 3*c->width }; // RGB stride
   sws_scale(ost->sws_ctx, (const uint8_t * const *)&RGB_data, inLinesize, 0, c->height, ost->frame->data, ost->frame->linesize);
   ost->frame->pts = ost->next_pts++;
   settings->frame_settings.zoom = settings->video_settings.zoom_func(ost->frame->pts);
   av_frame = ost->frame;

   if (av_frame == NULL) 
      return 1;
   av_init_packet(&pkt);
   ret = avcodec_send_frame(c, av_frame);
   if (ret == AVERROR_EOF)
      return 1;
   if (ret < 0) {
      fprintf(stderr, "Error on avcodec_send_frame: %s\n", av_err2str(ret));
      exit(1);
   }

   avcodec_receive_packet(c, &pkt);
   if (!ret)
      got_packet = 1;
   if (ret == AVERROR_EOF)
      return 1;
   if (ret == AVERROR(EAGAIN)) {
      fprintf(stderr, "Error on avcodec_receive_packet: %s\n", av_err2str(ret));
      return 0;
   }


   if (ret < 0) {
      fprintf(stderr, "Error encoding video frame: %s\n", av_err2str(ret));
      exit(1);
   }

   if (got_packet) {
      av_packet_rescale_ts(&pkt, c->time_base, ost->st->time_base);
      pkt.stream_index = ost->st->index;

      
      AVRational *time_base = &oc->streams[pkt.stream_index]->time_base;

      
      printf("Frame completed: pts:%s pts_time:%s dts:%s dts_time:%s\n",
             av_ts2str(pkt.pts), av_ts2timestr(pkt.pts, time_base),
             av_ts2str(pkt.dts), av_ts2timestr(pkt.dts, time_base));
       
      /* Write the compressed frame to the media file. */
      ret = av_interleaved_write_frame(oc, &pkt);
   } else {
      ret = 0;
   }
   if (ret == AVERROR_EOF)
      return 1;

   if (ret < 0) {
      fprintf(stderr, "Error while writing video frame: %s\n", av_err2str(ret));
      exit(1);
   }
   free_frame_data(&(frame_pointer->data));
   free(RGB_data);
   return (av_frame || got_packet) ? 0 : 1;
}


static void add_stream(OutputStream *ost, AVFormatContext *oc,
                       AVCodec **codec,
                       enum AVCodecID codec_id, const Settings *settings)
{
   AVCodecContext *c;
   /* find the encoder */
   *codec = avcodec_find_encoder(codec_id);
   if (!(*codec)) {
      fprintf(stderr, "Could not find encoder for '%s'\n", avcodec_get_name(codec_id));
      exit(1);
   }

   ost->st = avformat_new_stream(oc, NULL);
   if (!ost->st) {
      fprintf(stderr, "Could not allocate stream\n");
      exit(1);
   }
   ost->st->id = oc->nb_streams-1;
   c = avcodec_alloc_context3(*codec);
   if (!c) {
      fprintf(stderr, "Could not alloc an encoding context\n");
      exit(1);
   }
   ost->enc = c;
   c->codec_id = codec_id;
   c->bit_rate = 400000;
   /* Resolution must be a multiple of two. */
   c->width    = settings->mandelbrot_settings.f_res.x;
   c->height   = settings->mandelbrot_settings.f_res.x;
   /* timebase: This is the fundamental unit of time (in seconds) in terms
   * of which frame timestamps are represented. For fixed-fps content,
   * timebase should be 1/framerate and timestamp increments should be
   * identical to 1. */
   ost->st->time_base = (AVRational){ 1, settings->video_settings.frame_rate };
   c->time_base       = ost->st->time_base;

   c->gop_size      = 12; /* emit one intra frame every twelve frames at most */
   c->pix_fmt       = AV_PIX_FMT_YUV420P;
   c->max_b_frames = 2;

   /* Some formats want stream headers to be separate. */
   if (oc->oformat->flags & AVFMT_GLOBALHEADER)
      c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
}

static void close_stream(AVFormatContext *oc, OutputStream *ost) {
   avcodec_free_context(&ost->enc);
   av_frame_free(&ost->frame);
   av_frame_free(&ost->tmp_frame);
   sws_freeContext(ost->sws_ctx);
}

int render_video(const char *filename, Frame *frame_pointer, Settings *settings) {
   AVCodecContext *c= NULL;
   int ret;
   AVCodec *video_codec;
   OutputStream video_st = { 0 };
   AVOutputFormat *fmt;
   AVFormatContext *oc;
   AVDictionary *opt = NULL;

   av_dict_set(&opt, "x265-params", "qp=20", 0);
   av_dict_set(&opt, "preset", "medium", 0);
   av_dict_set(&opt, "tune", "zero-latency", 0);
   av_dict_set(&opt, "qmin", "0", 0);
   av_dict_set(&opt, "qmax", "69", 0);
   av_dict_set(&opt, "qdiff", "4", 0);

   av_register_all();
   avformat_alloc_output_context2(&oc, NULL, NULL, filename);
   if (!oc) {
      printf("Could not deduce output format from file extension: using MPEG.\n");
      avformat_alloc_output_context2(&oc, NULL, "mpeg", filename);
   }
   if (!oc) return 1;

   fmt = oc->oformat;
   fmt->video_codec = AV_CODEC_ID_HEVC;

   if (fmt->video_codec != AV_CODEC_ID_NONE) {
      add_stream(&video_st, oc, &video_codec, fmt->video_codec, settings);
   }

   c = video_st.enc;

   /* open it */
   avcodec_register(video_codec);
   ret = avcodec_open2(c, video_codec, &opt);
   av_dict_free(&opt);
   if (ret < 0) {
      fprintf(stderr, "Could not open codec: %s\n", av_err2str(ret));
      exit(1);
   }

    AVFrame *picture;
    picture = av_frame_alloc();
    if (!picture) {
      fprintf(stderr, "Could not allocate video frame\n");
      exit(1);
    }
    picture->format = c->pix_fmt;
    picture->width  = c->width;
    picture->height = c->height;

    /* allocate the buffers for the frame data */
    ret = av_frame_get_buffer(picture, 32);
    if (ret < 0) {
        fprintf(stderr, "Could not allocate frame data.\n");
        exit(1);
    }

   video_st.frame = picture;

    ret = avcodec_parameters_from_context(video_st.st->codecpar, c);
    if (ret < 0) {
        fprintf(stderr, "Could not copy the stream parameters\n");
        exit(1);
    }

   av_dump_format(oc, 0, filename, 1);

   if (!(fmt->flags & AVFMT_NOFILE)) {
      ret = avio_open(&oc->pb, filename, AVIO_FLAG_WRITE);
      if (ret < 0) {
         fprintf(stderr, "Could not open '%s': %s\n", filename,
         av_err2str(ret));
         return 1;
      }
    }

   /* Write the stream header, if any. */
   ret = avformat_write_header(oc, &opt);
   if (ret < 0) {
      fprintf(stderr, "Error occurred when opening output file: %s\n",
      av_err2str(ret));
      return 1;
   }

   int encode_video_ = 1;
   while (encode_video_) {
      encode_video_ = !write_video_frame(oc, &video_st, frame_pointer, settings);
   }


   av_write_trailer(oc);

   /* Close each codec. */
   close_stream(oc, &video_st);

   if (!(fmt->flags & AVFMT_NOFILE))
      avio_closep(&oc->pb);

   /* free the stream */
   avformat_free_context(oc);

   return 0;
}
