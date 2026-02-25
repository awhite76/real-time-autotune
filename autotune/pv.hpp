#ifndef PV_H
#define PV_H

#include "main.hpp"
#include <fftw3.h>

/* Functions to support phase vocoder implementation */

/* PV Params */
#define NUM_FRAMES PERIOD_FRAMES * CHANNELS // TOTAL number of samples per sampling period (PERIOD FRAMES is mono)
#define WINDOW_SIZE  NUM_FRAMES / 2
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

int settup_vocoder(float **time_buf, float **win, float **ifft_buf, float **omega, 
                    float **out, float **norm, int16_t **new_data, float **prev_phase, float **sum_phase, 
                    fftwf_complex **X, fftwf_complex **Y, 
                    float time_stretch, int* num_windows, int* Hs, 
                    int* out_L, fftwf_plan* p_r2c, fftwf_plan* p_c2r);

int phase_vocoder(int16_t* pcm, float *time_buf, float *win, float *ifft_buf, float* omega, 
                    float *out, float *norm, int16_t *new_data, float* prev_phase,float*  sum_phase, fftwf_complex *X, fftwf_complex *Y, 
                    int num_windows, int Hs, int out_L, fftwf_plan p_r2c, fftwf_plan p_c2r);





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