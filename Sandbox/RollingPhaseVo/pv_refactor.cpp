#include "pv_refactor.hpp"

int setup_vocoder(PhaseVocoder pv)
{
    if (!pv)
        return -1;

    // Start PV in known state
    memset(pv, 0, sizeof(*pv));

    // Ring sizes
    pv->input_length = (uint64_t)(8u * (uint64_t)WINDOW_SIZE);
    pv->stretched_length = (uint64_t)(64u * (uint64_t)WINDOW_SIZE);
    pv->norm_length = pv->stretched_length;

    // ---- FFT / time buffers ----
    pv->time_buf = (float *)fftwf_malloc(sizeof(float) * (size_t)WINDOW_SIZE);
    pv->ifft_buf = (float *)fftwf_malloc(sizeof(float) * (size_t)WINDOW_SIZE);
    pv->X = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (size_t)FREQ_BINS);
    pv->Y = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (size_t)FREQ_BINS);
    pv->win = (float *)fftwf_malloc(sizeof(float) * (size_t)WINDOW_SIZE);

    // ---- PV state ----
    pv->prev = (float *)calloc((size_t)FREQ_BINS, sizeof(float));
    pv->sum = (float *)calloc((size_t)FREQ_BINS, sizeof(float));
    pv->omega = (float *)malloc(sizeof(float) * (size_t)FREQ_BINS);

    // ---- Ring Buffers ----
    pv->stretched = (float *)calloc((size_t)pv->stretched_length, sizeof(float));
    pv->norm = (float *)calloc((size_t)pv->norm_length, sizeof(float));
    pv->input = (int16_t *)calloc((size_t)pv->input_length, sizeof(int16_t));

    if (!pv->time_buf || !pv->ifft_buf || !pv->X || !pv->Y || !pv->win ||
        !pv->prev || !pv->sum || !pv->omega || !pv->stretched || !pv->norm || !pv->input)
    {
        goto fail;
    }

    // Window
    _hann_window(pv->win, WINDOW_SIZE);

    // Plans
    pv->p_r2c = fftwf_plan_dft_r2c_1d(WINDOW_SIZE, pv->time_buf, pv->X, FFTW_ESTIMATE);
    pv->p_c2r = fftwf_plan_dft_c2r_1d(WINDOW_SIZE, pv->Y, pv->ifft_buf, FFTW_ESTIMATE);
    if (!pv->p_r2c || !pv->p_c2r)
    {
        goto fail;
    }

    // Expected bin center frequencies (rad/sample)
    for (int k = 0; k < FREQ_BINS; k++)
    {
        pv->omega[k] = 2.0f * (float)M_PI * (float)k / (float)WINDOW_SIZE;
    }

    pv->invN = 1.0f / (float)WINDOW_SIZE;

    pv->first_time = 1;

    pv->input_read = 0;
    pv->input_write = 0;

    pv->stretched_out_read = 0;
    pv->stretched_out_write = 0;

    return 0;

fail:
    // Plans
    if (pv->p_r2c)
        fftwf_destroy_plan(pv->p_r2c);
    if (pv->p_c2r)
        fftwf_destroy_plan(pv->p_c2r);

    // FFTW-allocated buffers
    if (pv->time_buf)
        fftwf_free(pv->time_buf);
    if (pv->ifft_buf)
        fftwf_free(pv->ifft_buf);
    if (pv->X)
        fftwf_free(pv->X);
    if (pv->Y)
        fftwf_free(pv->Y);
    if (pv->win)
        fftwf_free(pv->win);

    // malloc/calloc buffers
    free(pv->prev);
    free(pv->sum);
    free(pv->omega);
    free(pv->stretched);
    free(pv->norm);
    free(pv->input);

    // Leave pv in a clean state
    memset(pv, 0, sizeof(*pv));

    return -1;
}

void cleanup_vocoder(PhaseVocoder pv)
{
    if (!pv)
        return;

    // Destroy FFTW plans first
    if (pv->p_r2c)
    {
        fftwf_destroy_plan(pv->p_r2c);
        pv->p_r2c = NULL;
    }
    if (pv->p_c2r)
    {
        fftwf_destroy_plan(pv->p_c2r);
        pv->p_c2r = NULL;
    }

    // FFTW-allocated buffers
    if (pv->time_buf)
    {
        fftwf_free(pv->time_buf);
        pv->time_buf = NULL;
    }
    if (pv->ifft_buf)
    {
        fftwf_free(pv->ifft_buf);
        pv->ifft_buf = NULL;
    }
    if (pv->X)
    {
        fftwf_free(pv->X);
        pv->X = NULL;
    }
    if (pv->Y)
    {
        fftwf_free(pv->Y);
        pv->Y = NULL;
    }
    if (pv->win)
    {
        fftwf_free(pv->win);
        pv->win = NULL;
    }

    // malloc/calloc buffers
    free(pv->prev);
    pv->prev = NULL;
    free(pv->sum);
    pv->sum = NULL;
    free(pv->omega);
    pv->omega = NULL;
    free(pv->stretched);
    pv->stretched = NULL;
    free(pv->norm);
    pv->norm = NULL;
    free(pv->input);
    pv->input = NULL;

    // Reset scalar state
    memset(pv, 0, sizeof(*pv));
}

void phase_vocoder(PhaseVocoder pv)
{

    float *win = pv->win;
    float *time_buf = pv->time_buf;
    float *prev = pv->prev;
    float *sum = pv->sum;
    float *omega = pv->omega;
    float *ifft_buf = pv->ifft_buf;
    float *stretched = pv->stretched;
    float *norm = pv->norm;
    float invN = pv->invN;
    float time_stretch = pv->time_stretch;
    int16_t *input = pv->input;
    uint64_t stretched_length = pv->stretched_length;
    uint64_t stretched_out_write = pv->stretched_out_write;
    uint64_t input_read = pv->input_read;
    uint64_t input_length = pv->input_length;
    fftwf_plan p_r2c = pv->p_r2c;
    fftwf_plan p_c2r = pv->p_c2r;
    fftwf_complex *X = pv->X;
    fftwf_complex *Y = pv->Y;

    int Hs = (int)lroundf(ANALYSIS_HOP * time_stretch);

    /* Apply window to input */
    for (int n = 0; n < WINDOW_SIZE; n++)
    {
        int idx = (input_read + n) % (input_length);
        int16_t v = input[idx];
        float s = (float)v / 32768.0f;
        time_buf[n] = s * win[n];
    }

    /* Perform STFT */
    fftwf_execute(p_r2c);
    for (int k = 0; k < FREQ_BINS; k++)
    {
        const float re = X[k][0];
        float im = X[k][1];
        if (k == 0 || k == FREQ_BINS - 1)
        {
            im = 0;
        }
        const float mag = sqrtf(re * re + im * im);
        const float phase = atan2f(im, re);
        if (pv->first_time)
        {
            prev[k] = phase;
            sum[k] = phase;
        }
        else
        {
            // phase difference from expected
            float dphi = phase - prev[k] - omega[k] * (float)ANALYSIS_HOP;
            dphi = _princargf(dphi);
            // true frequency estimate (rad/sample)
            float true_freq = omega[k] + dphi / (float)ANALYSIS_HOP;
            // accumulate synthesis phase
            sum[k] += true_freq * (float)Hs;
            sum[k] = _princargf(sum[k]);
            prev[k] = phase;
        }

        Y[k][0] = mag * cosf(sum[k]);
        Y[k][1] = mag * sinf(sum[k]);
    }

    pv->first_time = 0;

    fftwf_execute(p_c2r);

    for (int n = 0; n < WINDOW_SIZE; n++)
    {

        int oidx = (stretched_out_write + n) % (stretched_length);
        float sample = ifft_buf[n] * invN;
        float wsample = sample * win[n];
        stretched[oidx] += wsample;
        norm[oidx] += win[n] * win[n];
    }

    (pv->stretched_out_write) = (stretched_out_write + Hs) % stretched_length;
    (pv->input_read) = (input_read + ANALYSIS_HOP) % input_length;
}

size_t pv_push_input(PhaseVocoder pv, const int16_t *buffer, size_t count)
{
    uint64_t input_read = pv->input_read;
    uint64_t input_write = pv->input_write;
    const uint64_t input_length = pv->input_length;

    size_t written = 0;
    while (written < count)
    {
        uint64_t free = _ring_free(input_read, input_write, input_length);
        if (free == 0)
            break;
        pv->input[input_write] = buffer[written];
        input_write = (input_write + 1) % input_length;
        written++;
    }

    pv->input_write = input_write;
    return written;
}

void pv_consume_output(PhaseVocoder pv, int16_t *out, size_t count)
{
    constexpr float EPS = 1e-12f;

    uint64_t stretched_out_read = pv->stretched_out_read;
    const uint64_t stretched_length = pv->stretched_length;

    for (size_t i = 0; i < count; i++)
    {
        const uint64_t idx = stretched_out_read;

        const float nrm = pv->norm[idx];
        float y = 0.0f;
        if (nrm > EPS)
        {
            y = pv->stretched[idx] / nrm;
        }

        // clear after read (critical!)
        pv->stretched[idx] = 0.0f;
        pv->norm[idx] = 0.0f;

        // float [-1,1] -> int16
        float s = std::max(-1.0f, std::min(1.0f, y));
        int32_t q = (int32_t)lrintf(s * 32767.0f);
        if (q < -32768)
            q = -32768;
        if (q > 32767)
            q = 32767;
        out[i] = (int16_t)q;

        stretched_out_read = (stretched_out_read + 1) % stretched_length;
    }

    pv->stretched_out_read = stretched_out_read;
}

size_t pv_process_ready(PhaseVocoder pv, int16_t *out, size_t out_cap)
{
    size_t produced = 0;
    uint64_t input_length = pv->input_length;

    while (_ring_used(pv->input_read, pv->input_write, input_length) >= (uint64_t)WINDOW_SIZE)
    {
        int Hs = (int)lroundf((float)ANALYSIS_HOP * pv->time_stretch);
        if (Hs <= 0)
            break;

        if (produced + (size_t)Hs > out_cap)
            break;

        // uint64_t out_free = _ring_free(pv->stretched_out_read,
        //                                pv->stretched_out_write,
        //                                pv->stretched_length);
        // if (out_free < (uint64_t)WINDOW_SIZE)
        // {
        //     // Not enough room to safely write the next iFFT frame.
        //     break;
        // }
        phase_vocoder(pv);
        pv_consume_output(pv, out + produced, (size_t)Hs);
        produced += (size_t)Hs;
    }

    return produced;
}
