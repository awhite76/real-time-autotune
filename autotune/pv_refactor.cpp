typedef struct PhaseVocoder_st {
    uint64_t in_window;   //
    uint64_t stretched_out_window;
    uint8_t first_time;
    float* win;
    float *prev_phase; 
    float *sum_phase; 
    float *omega;
    float *ifft_buf;
    float* time_buf;
    float time_stretch;
    float invN;
    fftwf_complex *X;
    fftwf_complex *Y;
    fftwf_plan p_r2c;
    fftwf_plan p_c2r;
    int16_t* input;
    boost::circular_buffer<float> norm;
    boost::circular_buffer<int16_t> out;
    boost::circular_buffer<float> stretched;
}* PhaseVocoder;

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

uint8_t phase_vocoder(PhaseVocoder pv) {

    uint64_t in_start = (pv->in_window) * ANALYSIS_HOP;
    /* Apply window to input */
    float s = 0.00;
    float* win = pv->win;
    float* time_buf = pv->time_buf;

    fftwf_plan p_r2c = pv->p_r2c;
    fftwf_plan p_c2r = pv->p_c2r;

    for (int n = 0; n < WINDOW_SIZE; n++) {
        int16_t v = input[n];
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

    int Hs = (int)lroundf(ANALYSIS_HOP * time_stretch);


    for (int k = 0; k < FREQ_BINS; k++) {

        const float re = X[k][0];
        const float im = X[k][1];
        const float mag = sqrtf(re*re + im*im);
        const float phase = atan2f(im, re);
        if (!pv->first_time) {
            prev[k] = phase;
            sum[k]  = phase;
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

    pv->first_time = 1;

    fftwf_execute(p_c2r);

    float invN = pv->invN;
    float* ifft_buf = pv->ifft_buf;
    uint64_t out_start = (pv->out_start) * Hs;

    for (int n = 0; n < WINDOW_SIZE; n++) {

        float sample = ifft_buf[n] * invN;
        float wsample = sample * win[n];

        stretched[out_start + n] += wsample;

        if (ch == 0) {
            // window-squared normalization for COLA robustness
            norm[out_start + n] += win[n] * win[n];
        }  
    }

    //TODO: Move pointers for input buffer to get the the next sample in the next call
    //TODO: Normalize the stretched ring buffer, then put it on the output.

    // Ensure all the pointers are moved appropriately
}