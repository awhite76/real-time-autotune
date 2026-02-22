#ifndef PV_H
#define PV_H

#include "main.hpp"
#include <fftw3.h>

/* Functions to support phase vocoder implementation */

/* Sampled data info */
#define NUM_CHANNELS 2
#define SAMPLING_RATE 48000
#define FRAME_BYTES 4
#define NUM_FRAMES BUFFER_FRAMES /* VERIFY THIS */

/* PV Params */
#define WINDOW_SIZE  480
#define ANALYSIS_HOP ((WINDOW_SIZE)/4)
#define FREQ_BINS (1 + WINDOW_SIZE/2) 

typedef struct {
    int N;              // FFT size / window length
    int H;              // hop size
    int num_frames;     // number of time frames
    float *window;      // length N
    float *time_buf;    // length N (input frame buffer)
    fftwf_complex *spec;// length num_frames * (N/2+1)
    fftwf_plan plan;    // r2c plan on time_buf -> temp_out
    fftwf_complex *temp_out; // length (N/2+1)
} STFT;



//float win[]
//float *time_buf, 
//float *ifft_buf, 
//fftwf_complex *X, 
//fftwf_complex *Y, 
//float *out, 
//float *norm,
//float *prev_phase, 
//float *sum_phase
//float *omega

#endif