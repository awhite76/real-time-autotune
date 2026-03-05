#ifndef RPV_H
#define RPV_H

#include <fftw3.h>
#include <stdint.h>
#include <algorithm>
#include <math.h>
#include <string.h>

typedef struct PhaseVocoder_st {
    int16_t *input;                  // input buffer (circular buffer)
    uint64_t input_read;             // Read SAMPLE for input
    uint64_t input_write;            // write SAMPLE for input
    uint64_t stretched_out_write;    // Write SAMPLE for stretched time output
    uint64_t stretched_out_read;     // Read SAMPLE for stretched time output
    uint64_t norm_length;            // length of normalization circular buffer
    uint64_t stretched_length;       // length of stretched time normalization buffer
    uint64_t input_length;           // Length of input ring buffer
    uint8_t first_time;              // Is this the first time that the PV is being called ?
    float *win;                      // Pointer to window
    float *prev;                     // Previous value at phase bin
    float *sum;                      // Accumulated value at phase bin
    float *omega;                    // Represents expected frequency at bins
    float *ifft_buf;                 // Output buffer of iFFT 
    float *time_buf;                 // Store windowed input
    float *norm;                     // Normalization buffer (Circular Buffer)
    float *stretched;                // Stretched output buffer (Circular Buffer)
    float time_stretch;              // ratio to time stretch (e.g. time_stretch of 2 is 2x time)
    float invN;                      // normalization factor
    fftwf_complex *X;                // Input buffer for FFT
    fftwf_complex *Y;                // Output buffer for results of iFFT
    fftwf_plan p_r2c;                // FFT plan
    fftwf_plan p_c2r;                // iFFT plan
} *PhaseVocoder;

static float _princargf(float x) {
    x = fmodf(x + (float)M_PI, 2.0f * (float)M_PI);
    if (x < 0) x += 2.0f * (float)M_PI;
    return x - (float)M_PI;
}

static void _hann_window(float *w, int N) {
    // periodic Hann: w[n] = 0.5 - 0.5*cos(2*pi*n/N)
    // (good for overlap-add with H=N/4 or H=N/2 depending)
    const float two_pi = 2.0f * (float)M_PI;
    for (int n = 0; n < N; n++) {
        w[n] = 0.5f - 0.5f * cosf(two_pi * (float)n / (float)N);
    }
}

// Used space
static inline uint64_t _ring_used(uint64_t read_pointer, uint64_t write_pointer, uint64_t ring_length) {
    return (write_pointer + ring_length - read_pointer) % ring_length;
}

// free space
static inline uint64_t _ring_free(uint64_t read_pointer, uint64_t write_pointer, uint64_t ring_length) {
    return (ring_length - 1) - _ring_used(read_pointer, write_pointer, ring_length);
}


int setup_vocoder(PhaseVocoder pv);
void cleanup_vocoder(PhaseVocoder pv);
size_t pv_process_ready(PhaseVocoder pv, int16_t* out, size_t out_cap);
void pv_consume_output(PhaseVocoder pv, int16_t* out, size_t count);
size_t pv_push_input(PhaseVocoder pv, const int16_t* buffer, size_t count); 

#define WINDOW_SIZE 1024
#define ANALYSIS_HOP 256
#define FREQ_BINS (1 + WINDOW_SIZE/2) 

#endif