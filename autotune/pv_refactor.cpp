typedef struct* PhaseVocoder_st {
    uint64_t in_window;
    uint64_t stretched_out_window;
    uint8_t first_time;
    float* win;
    float *prev; 
    float *sum; 
    float *omega;
    float *ifft_buf;
    float time_stretch;
    float invN;
    fftwf_complex *X;
    fftwf_complex *Y;
    fftwf_plan p_r2c;
    fftwf_plan p_c2r;
    boost::circular_buffer<int16_t> input(WINDOW_SIZE * 5);
    boost::circular_buffer<int16_t> stretched(WINDOW_SIZE * 5);
    boost::circular_buffer<int16_t> norm(WINDOW_SIZE * 5);
} PhaseVocoder;

float _princargf(float x) {
    x = fmodf(x + (float)M_PI, 2.0f * (float)M_PI);
    if (x < 0) x += 2.0f * (float)M_PI;
    return x - (float)M_PI;
}


//TODO
uint8_t settup_vocoder() {
    return -1;
}

//TODO
uint8_t cleanup_vocoder() {
    return -1;
}

uint8_t phase_vocoder(input_windowPhaseVocoder pv) {

    /* Want to input a means of figuring out starting point in ring buffer....
            then that would also work for the norm one too */

    uint64_t in_start = (pv->in_start) * ANALYSIS_HOP;
    /* Apply window to input */
    float s = 0.00;
    float* win = pv->win;
    for (int n = 0; n < WINDOW_SIZE; n++) {
        int idx = in_start + n;
        int16_t v = input[idx];
        s = (float)v / 32768.0f;
        time_buf[n] = s * win[n];
    }

    /* Perform STFT */
    fftwf_execute(p_r2c);

    // Make these typedef's?
    fftwf_complex* X = pv->X;
    fftwf_complex* Y = pv->Y;
    float* prev  = pv->prev;
    float* sum   = pv->sum;
    float* omega = pv->omega;

    for (int k = 0; k < FREQ_BINS; k++) {
        const float re = X[k][0];
        const float im = X[k][1];
        const float mag = sqrtf(re*re + im*im);
        const float phase = atan2f(im, re);
        if (!pv->first_time) {
            prev[k] = phase;
            sum[k]  = phase;
            pv->first_time = 1;
        } else {
            // phase difference from expected
            float dphi = phase - prev[k] - omega[k] * (float)ANALYSIS_HOP;
            dphi = _princargf(dphi)
            // true frequency estimate (rad/sample)
            float true_freq = omega[k] + dphi / (float)ANALYSIS_HOP
            // accumulate synthesis phase
            sum[k] += true_freq * (float)Hs;
            prev[k] = phase;
        
            Y[k][0] = mag * cosf(sum[k]);
            Y[k][1] = mag * sinf(sum[k]);
        }
    }

    fftwf_execute(p_c2r);

    float invN = pv->invN;
    float* ifft_buf = pv->ifft_buf;

    int Hs = (int)lroundf(ANALYSIS_HOP * time_stretch);

    uint64_t out_start = (pv->out_start) * Hs;

    static uint8_t ch = 0;
    for (int n = 0; n < WINDOW_SIZE; n++) {
        int oidx = out_start + n;

        float sample = ifft_buf[n] * invN;
        float wsample = sample * win[n];

        out[oidx + n] += wsample;

        if (ch == 0) {
            // window-squared normalization for COLA robustness
            norm[oidx] += win[n] * win[n];
            ch = ch - 1;
        }  
    }

    for (int n = 0; n < WINDOW_SIZE; n++) {
        int oidx = out_start + n;
        float g = norm[oidx];
        if (g < 1e-12f) g = 1.0f;
        float invg = 1.0f / g;
        out[oidx] *= invg;
    }

    for (int n = 0; n < WINDOW_SIZE; n++) {
            int oidx = out_start + n;
            float v = out[oidx];
            // simple clip
            if (v > 1.0f) v = 1.0f;
            if (v < -1.0f) v = -1.0f;
            int32_t q = (int32_t)lroundf(v * 32767.0f);
            if (q > 32767) q = 32767;
            if (q < -32768) q = -32768;
            new_data[(size_t)n * NUM_CHANNELS + (size_t)ch] = (int16_t)q;

    }
}