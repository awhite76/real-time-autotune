int phase_vocoder(int16_t* pcm, float *time_buf, float *win, float *ifft_buf, float* omega, 
                  float *out, float *norm, int16_t *new_data, float* prev_phase, float*  sum_phase, fftwf_complex *X, fftwf_complex *Y, 
                  float time_stretch, int *out_L, int num_windows, fftwf_plan p_r2c, fftwf_plan p_c2r) {

    
    if(time_stretch > MAX_TIME_STRETCH) {
        time_stretch = MAX_TIME_STRETCH;
    }
    int Hs = (int)lroundf(ANALYSIS_HOP * time_stretch);
    *out_L = (num_windows - 1) * (Hs) + WINDOW_SIZE; 

    // --------------------------
    // Perform PhaseVo Algorithm on each channel
    // --------------------------
    for (int ch = 0; ch < NUM_CHANNELS; ch++) {
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
                    int16_t v = pcm[(size_t)idx * NUM_CHANNELS + (size_t)ch];
                    s = (float)v / 32768.0f;
                }
                time_buf[n] = s * win[n];
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
                    // phase difference from expected
                    float dphi = phase - prev[k] - omega[k] * (float)ANALYSIS_HOP;
                    dphi = princargf(dphi);

                    // true frequency estimate (rad/sample)
                    float true_freq = omega[k] + dphi / (float)ANALYSIS_HOP;

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
            const float invN = 1.0f / (float)WINDOW_SIZE;

            // overlap-add (window again) + norm accumulation once (for ch==0)
            for (int n = 0; n < WINDOW_SIZE; n++) {
                int oidx = out_start + n;
                if (oidx >= *out_L) break;

                float sample = ifft_buf[n] * invN;
                float wsample = sample * win[n];

                out[(size_t)oidx * (size_t)NUM_CHANNELS + (size_t)ch] += wsample;

                if (ch == 0) {
                    // window-squared normalization for COLA robustness
                    norm[oidx] += win[n] * win[n];
                }
            }
        }
    }

    // Normalize overlap-add
    for (int n = 0; n < *out_L; n++) {
        float g = norm[n];
        if (g < 1e-12f) g = 1.0f;
        float invg = 1.0f / g;
        for (int ch = 0; ch < NUM_CHANNELS; ch++) {
            out[(size_t)n * NUM_CHANNELS + (size_t)ch] *= invg;
        }
    }
    
    for (int n = 0; n < *out_L; n++) {
        for (int ch = 0; ch < NUM_CHANNELS; ch++) {
            float v = out[(size_t)n * NUM_CHANNELS + (size_t)ch];
            // simple clip
            if (v > 1.0f) v = 1.0f;
            if (v < -1.0f) v = -1.0f;
            int32_t q = (int32_t)lroundf(v * 32767.0f);
            if (q > 32767) q = 32767;
            if (q < -32768) q = -32768;
            new_data[(size_t)n * NUM_CHANNELS + (size_t)ch] = (int16_t)q;
        }
    }

    return 0;
}
