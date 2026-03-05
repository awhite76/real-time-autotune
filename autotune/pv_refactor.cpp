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

void phase_vocoder(PhaseVocoder pv) {

    float* win = pv->win;
    float* time_buf = pv->time_buf;
    float* prev  = pv->prev;
    float* sum   = pv->sum;
    float* omega = pv->omega;
    float* ifft_buf = pv->ifft_buf;
    float* stretched = pv->stretched;
    float* norm = pv->norm;
    float invN = pv->invN;
    float time_stretch = pv->time_stretch;
    int16_t* input = pv->input;
    uint64_t stretched_length = pv->stretched_length;
    uint64_t stretched_out_write = pv->stretched_out_write;
    uint64_t input_read = pv->input_read;
    uint64_t input_length = pv->input_length;
    fftwf_plan p_r2c = pv->p_r2c;
    fftwf_plan p_c2r = pv->p_c2r;
    fftwf_complex* X = pv->X;
    fftwf_complex* Y = pv->Y;

    int Hs = (int)lroundf(ANALYSIS_HOP * time_stretch);

    /* Apply window to input */
    for (int n = 0; n < WINDOW_SIZE; n++) {
        int idx = (input_read + n) % (input_length);
        int16_t v = input[idx];
        float s = (float)v / 32768.0f;
        time_buf[n] = s * win[n];
    }

    /* Perform STFT */
    fftwf_execute(p_r2c);
    for (int k = 0; k < FREQ_BINS; k++) {
        const float re = X[k][0];
        const float im = X[k][1];
        const float mag = sqrtf(re*re + im*im);
        const float phase = atan2f(im, re);
        if (pv->first_time) {
            prev[k] = phase;
            sum[k]  = phase;
        } else {
            // phase difference from expected
            float dphi = phase - prev[k] - omega[k] * (float)ANALYSIS_HOP;
            dphi = _princargf(dphi);
            // true frequency estimate (rad/sample)
            float true_freq = omega[k] + dphi / (float)ANALYSIS_HOP;
            // accumulate synthesis phase
            sum[k] += true_freq * (float)Hs;
            prev[k] = phase;
        
        }

        Y[k][0] = mag * cosf(sum[k]);
        Y[k][1] = mag * sinf(sum[k]);
    }

    pv->first_time = 0;

    fftwf_execute(p_c2r);

    for (int n = 0; n < WINDOW_SIZE; n++) {
        
        int oidx = (stretched_out_write + n) % (stretched_length);
        float sample = ifft_buf[n] * invN;
        float wsample = sample * win[n];
        stretched[oidx] += wsample;
        norm[oidx] += win[n] * win[n]; 
    }

    (pv->stretched_out_write) = (stretched_out_write + Hs) % stretched_length;
    (pv->input_read) = (input_read+ANALYSIS_HOP) % input_length;
}