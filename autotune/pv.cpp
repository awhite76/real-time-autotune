#include "pv.hpp"
#include <math.h>


float princargf(float x) {
    // wrap to (-pi, pi]
    x = fmodf(x + (float)M_PI, 2.0f * (float)M_PI);
    if (x < 0) x += 2.0f * (float)M_PI;
    return x - (float)M_PI;
}

float fast_atan2f(float y, float x) {
    return atan2f(y, x);
}

void hann_window(float *w, int N) {
    // periodic Hann: w[n] = 0.5 - 0.5*cos(2*pi*n/N)
    const float two_pi = 2.0f * (float)M_PI;
    for (int n = 0; n < N; n++) {
        w[n] = 0.5f - 0.5f * cosf(two_pi * (float)n / (float)N);
    }
}

int settup_vocoder(float **time_buf, float **win, float **ifft_buf, float **omega, 
                    float **out, float **norm, int16_t **new_data, float **prev_phase, float **sum_phase, 
                    fftwf_complex **X, fftwf_complex **Y, int* num_windows, int* Hs, int* out_L, fftwf_plan* p_r2c, fftwf_plan* p_c2r)

{

    *time_buf = (float*) fftwf_malloc(sizeof(float) * (size_t)WINDOW_SIZE);
    *ifft_buf = (float*) fftwf_malloc(sizeof(float) * (size_t)WINDOW_SIZE);
    *X = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)FREQ_BINS);
    *Y = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)FREQ_BINS);
    if (!*time_buf || !*ifft_buf || !*X || !*Y) {
        if (*time_buf) fftwf_free(*time_buf);
        if (*ifft_buf) fftwf_free(*ifft_buf);
        if (*X) fftwf_free(*X);
        if (*Y) fftwf_free(*Y);
        return -1;
    }

    *Hs = (int)lroundf(ANALYSIS_HOP * MAX_TIME_STRETCH);
    *num_windows = 1 + (int)ceilf((float)((NUM_FRAMES) - WINDOW_SIZE) / (float)ANALYSIS_HOP);
    *out_L = (*num_windows - 1) * (*Hs) + WINDOW_SIZE; 
    *win = (float*)fftwf_malloc(sizeof(float) * (size_t)WINDOW_SIZE);
    if (!(*win)) return -1;
    hann_window(*win, WINDOW_SIZE);

    *p_r2c = fftwf_plan_dft_r2c_1d(WINDOW_SIZE, *time_buf, *X, FFTW_ESTIMATE);

    *p_c2r = fftwf_plan_dft_c2r_1d(WINDOW_SIZE, *Y, *ifft_buf, FFTW_ESTIMATE);

    if (!(*p_r2c) || !(*p_c2r)) {
        if (*p_r2c) fftwf_destroy_plan(*p_r2c);
        if (*p_c2r) fftwf_destroy_plan(*p_c2r);
        return -1;
    }

    *out = (float*)calloc((size_t)(*out_L) * (size_t)CHANNELS, sizeof(float));
    *norm = (float*)calloc((size_t)(*out_L), sizeof(float)); // same for all channels
    if (!(*out) || !(*norm)) {
        free(*out); free(*norm);
        return -1;
    }

    // Calculate normalization buffer for IFFT
    for(int n = 0; n < WINDOW_SIZE; n++) {
        norm[n] = win[n]*win[n];    
    }

    *prev_phase = (float*)calloc((size_t)CHANNELS * (size_t)FREQ_BINS, sizeof(float));
    *sum_phase  = (float*)calloc((size_t)CHANNELS * (size_t)FREQ_BINS, sizeof(float));
    if (!*prev_phase || !*sum_phase) {
        free(*prev_phase); free(*sum_phase);
        free(*out); free(*norm);
        fftwf_destroy_plan(*p_r2c); fftwf_destroy_plan(*p_c2r);
        fftwf_free(*win); fftwf_free(*time_buf); fftwf_free(*ifft_buf); fftwf_free(*X); fftwf_free(*Y);
        return -1;
    }

    *omega = (float*)malloc(sizeof(float) * (size_t)FREQ_BINS);
    if (!(*omega)) {
        free(*prev_phase); free(*sum_phase);
        free(*out); free(*norm);
        fftwf_destroy_plan(*p_r2c); fftwf_destroy_plan(*p_c2r);
        fftwf_free(*win); fftwf_free(*time_buf); fftwf_free(*ifft_buf); fftwf_free(*X); fftwf_free(*Y);
        return -1;
    }

    for (int k = 0; k < FREQ_BINS; k++) {
        (*omega)[k] = 2.0f * (float)M_PI * (float)k / (float)WINDOW_SIZE; // radians/sample
    }

    *new_data = (int16_t*)malloc((*out_L) * CHANNELS * sizeof(int16_t));
    if (!*new_data) {
        free(*omega);
        free(*prev_phase); free(*sum_phase);
        free(*out); free(*norm);
        fftwf_destroy_plan(*p_r2c); fftwf_destroy_plan(*p_c2r);
        fftwf_free(*win); fftwf_free(*time_buf); fftwf_free(*ifft_buf); fftwf_free(*X); fftwf_free(*Y);
        return -1;
    }

    return 0;
}

// void cleanup_vocoder() {
//     //TODO:
// }



int phase_vocoder(int16_t* pcm, float *time_buf, float *win, float *ifft_buf, float* omega, 
                  float *out, float *norm, int16_t *new_data, float* prev_phase, float*  sum_phase, fftwf_complex *X, fftwf_complex *Y, 
                  int num_windows, float time_stretch, fftwf_plan p_r2c, fftwf_plan p_c2r) {

    // Calculate analysis hop 

    if(time_stretch > MAX_TIME_STRETCH) {
        time_stretch =  MAX_TIME_STRETCH;
    }

    int Hs = (int)lroundf(ANALYSIS_HOP * time_stretch);

    int out_L = (num_windows - 1) * (Hs) + WINDOW_SIZE;


    // --------------------------
    // Perform PhaseVo Algorithm on each channel
    // --------------------------
    for (int ch = 0; ch < CHANNELS; ch++) {
        float *prev = prev_phase + (size_t)ch * FREQ_BINS;
        float *sum  = sum_phase  + (size_t)ch * FREQ_BINS;

        // init phases from first frame
        // (we’ll do it naturally on f=0)
        for (int f = 0; f < num_windows; f++) {
            const int in_start  = f * ANALYSIS_HOP;
            const int out_start = f * Hs;

            // Analysis frame: windowed input into time_buf
            for (int n = 0; n < WINDOW_SIZE; n++) {
                int idx = in_start + n;
                float s = 0.0f;
                if (idx < NUM_FRAMES) {
                    // interleaved: frame idx, channel ch
                    int16_t v = pcm[(size_t)idx * CHANNELS + (size_t)ch];
                    s = (float)v / 32768.0f;

                    //SIMD here ?
                }
                time_buf[n] = s * win[n];

                //SIMD here ?
            }

            fftwf_execute(p_r2c);

            // Build Y with phase propagation
            for (int k = 0; k < FREQ_BINS; k++) {
                const float re = X[k][0];
                const float im = X[k][1];

                const float mag = sqrtf(re*re + im*im);
                const float phase = fast_atan2f(im, re);

                if (f == 0) {
                    prev[k] = phase;
                    sum[k]  = phase;
                } else {

                    // SIMD here ?
                    // phase difference from expected
                    float dphi = phase - prev[k] - omega[k] * (float)ANALYSIS_HOP;
                    dphi = princargf(dphi);

                    // true frequency estimate (rad/sample)
                    float true_freq = omega[k] + dphi / (float)ANALYSIS_HOP;

                    // accumulate synthesis phase
                    sum[k] += true_freq * (float)Hs;
                    prev[k] = phase;
                }

                //SIMD here ?
                Y[k][0] = mag * cosf(sum[k]);
                Y[k][1] = mag * sinf(sum[k]);
            }

            // ISTFT
            fftwf_execute(p_c2r);

            // FFTW's c2r is unnormalized: divide by N
            const float invN = 1.0f / (float)WINDOW_SIZE;

            // overlap-add (window again) + norm accumulation once (for ch==0)
            for (int n = 0; n < WINDOW_SIZE; n++) {

                // SIMD here ?
                int oidx = out_start + n;
                if (oidx >= out_L) break;

                float sample = ifft_buf[n] * invN;
                float wsample = sample * win[n];

                out[(size_t)oidx * (size_t)CHANNELS + (size_t)ch] += wsample;

                if (ch == 0) {
                    // window-squared normalization for COLA robustness
                    norm[oidx] += win[n] * win[n];
                }
            }
        }
    }

    // Normalize overlap-add
    for (int n = 0; n < out_L; n++) {

        //SIMD here ?
        float g = norm[n];
        if (g < 1e-12f) g = 1.0f;
        float invg = 1.0f / g;
        for (int ch = 0; ch < CHANNELS; ch++) {
            out[(size_t)n * CHANNELS + (size_t)ch] *= invg;
        }
    }
    
    for (int n = 0; n < out_L; n++) {
        for (int ch = 0; ch < CHANNELS; ch++) {
            float v = out[(size_t)n * CHANNELS + (size_t)ch];
            // simple clip
            if (v > 1.0f) v = 1.0f;
            if (v < -1.0f) v = -1.0f;
            int32_t q = (int32_t)lroundf(v * 32767.0f);
            if (q > 32767) q = 32767;
            if (q < -32768) q = -32768;
            new_data[(size_t)n * CHANNELS + (size_t)ch] = (int16_t)q;
        }
    }

    return 0;
}
