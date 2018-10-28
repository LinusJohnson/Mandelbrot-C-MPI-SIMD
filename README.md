This program renders images and a fractal zoom animation of the Mandelbrot set.

Compile using make and then run using the flags:

-n to run in non-SIMD mode (SIMD mode on by default)  
-v n to render a video with n frames  
-f n to specify the frame rate (default 30)  
-r n to set the resolution to n x n (default 2048 x 2048)  
-t n to run in performance test mode with n tests (must have been compiled with TIME_FRAMES set to 1 in definitions.h)
