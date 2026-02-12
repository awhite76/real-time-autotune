#include "pv.hpp"


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
    // (good for overlap-add with H=N/4 or H=N/2 depending)
    const float two_pi = 2.0f * (float)M_PI;
    for (int n = 0; n < N; n++) {
        w[n] = 0.5f - 0.5f * cosf(two_pi * (float)n / (float)N);
    }
}

int phase_vocoder(int16_t* buffer, float time_stretch) {

    /* Get data from wav file */
    const int C = (int)w->fmt.num_channels;
    const int sr = (int)w->fmt.sample_rate;
    const uint32_t frame_bytes = w->fmt.block_align; // C * 2
    const int L = (int)(w->data_bytes / frame_bytes); // number of frames (samples per channel)
    const int16_t *pcm = (const int16_t*)w->data;

    /* Phase Vocoder parameters */
    const int N  = 1024;   // FFT / window size
    const int Ha = 256;    // analysis hop (75% overlap)
    const int Hs = (int)lroundf((float)Ha * time_stretch); // synthesis hop

    const int K = N/2 + 1;
    //const int num_windows = stft_num_frames(L, N, Ha);
    const int num_windows = 1 + (int)ceilf((float)(L - N) / (float)Ha); // number of windows for the STFT
    const int out_L = (num_windows - 1) * Hs + N; // output samples per channel

    // window
    //float *win = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    //if (!win) return -1;
    //hann_window(win, N);

    // FFTW buffers/plans
    //float *time_buf = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    //float *ifft_buf = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    //fftwf_complex *X = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)K);
    //fftwf_complex *Y = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)K);
    //if (!time_buf || !ifft_buf || !X || !Y) {
    //    fftwf_free(win);
    //    if (time_buf) fftwf_free(time_buf);
    //    if (ifft_buf) fftwf_free(ifft_buf);
    //    if (X) fftwf_free(X);
    //    if (Y) fftwf_free(Y);
    //    return -1;
    //}

    fftwf_plan p_r2c = fftwf_plan_dft_r2c_1d(N, time_buf, X, FFTW_MEASURE);
    fftwf_plan p_c2r = fftwf_plan_dft_c2r_1d(N, Y, ifft_buf, FFTW_MEASURE);
    if (!p_r2c || !p_c2r) {
        if (p_r2c) fftwf_destroy_plan(p_r2c);
        if (p_c2r) fftwf_destroy_plan(p_c2r);
        return -1;
    }

    // Output buffers per channel
    //float *out = (float*)calloc((size_t)out_L * (size_t)C, sizeof(float));
    //float *norm = (float*)calloc((size_t)out_L, sizeof(float)); // same for all channels
    //if (!out || !norm) {
    //    free(out); free(norm);
    //    fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
    //    fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
    //    return -1;
    //}

    // Phase vocoder state per channel per bin
    //float *prev_phase = (float*)calloc((size_t)C * (size_t)K, sizeof(float));
    //float *sum_phase  = (float*)calloc((size_t)C * (size_t)K, sizeof(float));
    //if (!prev_phase || !sum_phase) {
    //    free(prev_phase); free(sum_phase);
    //    free(out); free(norm);
    //    fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
    //    fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
    //    return -1;
    //}

    // Precompute expected phase advance per bin for analysis hop
    //float *omega = (float*)malloc(sizeof(float) * (size_t)K);
    //if (!omega) {
    //    free(prev_phase); free(sum_phase);
    //    free(out); free(norm);
    //    fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
    //    fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
    //    return -1;
    //}
    for (int k = 0; k < K; k++) {
        omega[k] = 2.0f * (float)M_PI * (float)k / (float)N; // radians/sample
    }

    // --------------------------
    // Process each channel independently
    // --------------------------
    for (int ch = 0; ch < C; ch++) {
        float *prev = prev_phase + (size_t)ch * (size_t)K;
        float *sum  = sum_phase  + (size_t)ch * (size_t)K;

        // init phases from first frame
        // (we’ll do it naturally on f=0)
        for (int f = 0; f < num_windows; f++) {
            const int in_start  = f * Ha;
            const int out_start = f * Hs;

            // Analysis frame: windowed input into time_buf
            for (int n = 0; n < N; n++) {
                int idx = in_start + n;
                float s = 0.0f;
                if (idx < L) {
                    // interleaved: frame idx, channel ch
                    int16_t v = pcm[(size_t)idx * (size_t)C + (size_t)ch];
                    s = (float)v / 32768.0f;
                }
                time_buf[n] = s * win[n];
            }

            fftwf_execute(p_r2c);

            // Build Y with phase propagation
            for (int k = 0; k < K; k++) {
                const float re = X[k][0];
                const float im = X[k][1];

                const float mag = sqrtf(re*re + im*im);
                const float phase = fast_atan2f(im, re);

                if (f == 0) {
                    prev[k] = phase;
                    sum[k]  = phase;
                } else {
                    // phase difference from expected
                    float dphi = phase - prev[k] - omega[k] * (float)Ha;
                    dphi = princargf(dphi);

                    // true frequency estimate (rad/sample)
                    float true_freq = omega[k] + dphi / (float)Ha;

                    // accumulate synthesis phase
                    sum[k] += true_freq * (float)Hs;
                    prev[k] = phase;
                }

                Y[k][0] = mag * cosf(sum[k]);
                Y[k][1] = mag * sinf(sum[k]);
            }

            // ISTFT
            fftwf_execute(p_c2r);

            // FFTW's c2r is unnormalized: divide by N
            const float invN = 1.0f / (float)N;

            // overlap-add (window again) + norm accumulation once (for ch==0)
            for (int n = 0; n < N; n++) {
                int oidx = out_start + n;
                if (oidx >= out_L) break;

                float sample = ifft_buf[n] * invN;
                float wsample = sample * win[n];

                out[(size_t)oidx * (size_t)C + (size_t)ch] += wsample;

                if (ch == 0) {
                    // window-squared normalization for COLA robustness
                    norm[oidx] += win[n] * win[n];
                }
            }
        }
    }

    // Normalize overlap-add
    for (int n = 0; n < out_L; n++) {
        float g = norm[n];
        if (g < 1e-12f) g = 1.0f;
        float invg = 1.0f / g;
        for (int ch = 0; ch < C; ch++) {
            out[(size_t)n * (size_t)C + (size_t)ch] *= invg;
        }
    }

    // Convert back to int16 PCM interleaved
    const uint32_t out_data_bytes = (uint32_t)((size_t)out_L * (size_t)C * sizeof(int16_t));
    uint8_t *new_data = (uint8_t*)malloc(out_data_bytes);
    if (!new_data) {
        // cleanup and leave wav unchanged
        free(omega);
        free(prev_phase); free(sum_phase);
        free(out); free(norm);
        fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
        fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
        return -1;
    }

    int16_t *out_pcm = (int16_t*)new_data;
    for (int n = 0; n < out_L; n++) {
        for (int ch = 0; ch < C; ch++) {
            float v = out[(size_t)n * (size_t)C + (size_t)ch];
            // simple clip
            if (v > 1.0f) v = 1.0f;
            if (v < -1.0f) v = -1.0f;
            int32_t q = (int32_t)lroundf(v * 32767.0f);
            if (q > 32767) q = 32767;
            if (q < -32768) q = -32768;
            out_pcm[(size_t)n * (size_t)C + (size_t)ch] = (int16_t)q;
        }
    }

    // Replace wav data
    free(w->data);
    w->data = new_data;
    w->data_bytes = out_data_bytes;

    // byte_rate / block_align / sample_rate / channels remain correct for PCM
    w->fmt.byte_rate = w->fmt.sample_rate * w->fmt.block_align;

    // cleanup
    free(omega);
    free(prev_phase); free(sum_phase);
    free(out); free(norm);
    fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
    fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);

    fprintf(stderr, "phase_vocoder: sr=%d, ch=%d, in_L=%d, out_L=%d, stretch=%.3f\n",
            sr, C, L, out_L, time_stretch);

    return 0;
}