// PV.c
// Build: gcc -O2 -Wall PV.c -o pv
// Run:   ./pv in.wav out.wav

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>

/* Initial settup */

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

#pragma pack(push, 1)
typedef struct {
    char     riff_id[4];   // "RIFF"
    uint32_t riff_size;    // file_size - 8
    char     wave_id[4];   // "WAVE"
} RiffHeader;

typedef struct {
    char     id[4];
    uint32_t size;
} ChunkHeader;

typedef struct {
    uint16_t audio_format;    // 1 = PCM
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bits_per_sample;
    // may be followed by extra fmt bytes if chunk size > 16
} FmtChunk;
#pragma pack(pop)

typedef struct {
    FmtChunk fmt;
    uint8_t *data;          // raw PCM bytes
    uint32_t data_bytes;    // bytes in data chunk
} WavFile;

static int read_exact(FILE *f, void *buf, size_t n) {
    return fread(buf, 1, n, f) == n ? 0 : -1;
}

static int write_exact(FILE *f, const void *buf, size_t n) {
    return fwrite(buf, 1, n, f) == n ? 0 : -1;
}

static void free_wav(WavFile *w) {
    free(w->data);
    memset(w, 0, sizeof(*w));
}

/* For SFTF */

static void hann_window(float *w, int N) {
    // periodic Hann: w[n] = 0.5 - 0.5*cos(2*pi*n/N)
    // (good for overlap-add with H=N/4 or H=N/2 depending)
    const float two_pi = 2.0f * (float)M_PI;
    for (int n = 0; n < N; n++) {
        w[n] = 0.5f - 0.5f * cosf(two_pi * (float)n / (float)N);
    }
}

static int stft_init(STFT *s, int N, int H, int num_frames) {
    memset(s, 0, sizeof(*s));
    s->N = N;
    s->H = H;
    s->num_frames = num_frames;

    s->window = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    s->time_buf = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    s->temp_out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)(N/2 + 1));
    s->spec = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)num_frames * (size_t)(N/2 + 1));

    if (!s->window || !s->time_buf || !s->temp_out || !s->spec) return -1;

    hann_window(s->window, N);

    // FFTW plan: real input -> complex output
    s->plan = fftwf_plan_dft_r2c_1d(N, s->time_buf, s->temp_out, FFTW_MEASURE);
    if (!s->plan) return -1;

    return 0;
}

static void stft_free(STFT *s) {
    if (s->plan) fftwf_destroy_plan(s->plan);
    if (s->window) fftwf_free(s->window);
    if (s->time_buf) fftwf_free(s->time_buf);
    if (s->temp_out) fftwf_free(s->temp_out);
    if (s->spec) fftwf_free(s->spec);
    memset(s, 0, sizeof(*s));
}


/* Make test sine wave */
static void make_test_sine(float *x, int L, float sr, float freq) {
    const float w = 2.0f * (float)M_PI * freq / sr;
    for (int n = 0; n < L; n++) x[n] = sinf(w * (float)n);
}

// Compute how many frames you’ll produce for a signal of length L
// This version includes a final zero-padded frame if needed.
static int stft_num_frames(int L, int N, int H) {
    if (L <= 0) return 0;
    if (L <= N) return 1;
    int frames = 1 + (int)ceilf((float)(L - N) / (float)H);
    return frames;
}

// STFT: input x[0..L-1] -> s->spec
// spec is stored frame-major: spec[frame*(N/2+1) + bin]
static void stft_compute(STFT *s, const float *x, int L) {
    const int N = s->N;
    const int H = s->H;
    const int K = N/2 + 1;

    for (int f = 0; f < s->num_frames; f++) {
        int start = f * H;

        // gather frame with zero-padding
        for (int n = 0; n < N; n++) {
            int idx = start + n;
            float sample = (idx < L) ? x[idx] : 0.0f;
            s->time_buf[n] = sample * s->window[n];
        }

        // FFT
        fftwf_execute(s->plan);

        // copy to output array
        fftwf_complex *dst = s->spec + (size_t)f * (size_t)K;
        memcpy(dst, s->temp_out, sizeof(fftwf_complex) * (size_t)K);
    }
}

static int load_wav_pcm(const char *path, WavFile *out) {
    memset(out, 0, sizeof(*out));

    FILE *f = fopen(path, "rb");
    if (!f) { perror("fopen"); return -1; }

    RiffHeader rh;
    if (read_exact(f, &rh, sizeof(rh)) != 0) { fclose(f); return -1; }

    if (memcmp(rh.riff_id, "RIFF", 4) != 0 || memcmp(rh.wave_id, "WAVE", 4) != 0) {
        fprintf(stderr, "Not a RIFF/WAVE file\n");
        fclose(f);
        return -1;
    }

    int have_fmt = 0, have_data = 0;
    FmtChunk fmt = {0};
    uint8_t *data = NULL;
    uint32_t data_bytes = 0;

    while (!have_data) {
        ChunkHeader ch;
        if (read_exact(f, &ch, sizeof(ch)) != 0) break;

        if (memcmp(ch.id, "fmt ", 4) == 0) {
            if (ch.size < 16) { fprintf(stderr, "fmt chunk too small\n"); fclose(f); return -1; }
            if (read_exact(f, &fmt, 16) != 0) { fclose(f); return -1; }
            if (ch.size > 16) fseek(f, (long)(ch.size - 16), SEEK_CUR);
            have_fmt = 1;
        } else if (memcmp(ch.id, "data", 4) == 0) {
            if (!have_fmt) { fprintf(stderr, "data before fmt (unsupported here)\n"); fclose(f); return -1; }
            data = (uint8_t*)malloc(ch.size);
            if (!data) { fclose(f); return -1; }
            if (read_exact(f, data, ch.size) != 0) { free(data); fclose(f); return -1; }
            data_bytes = ch.size;
            have_data = 1;
        } else {
            // skip unknown chunk
            fseek(f, (long)ch.size, SEEK_CUR);
        }

        // chunks are word-aligned: pad byte if odd size
        if (ch.size & 1) fseek(f, 1, SEEK_CUR);
    }

    fclose(f);

    if (!have_fmt || !have_data) {
        fprintf(stderr, "Missing fmt or data chunk\n");
        free(data);
        return -1;
    }
    if (fmt.audio_format != 1) {
        fprintf(stderr, "Unsupported WAV format %u (supports PCM=1 only)\n", fmt.audio_format);
        free(data);
        return -1;
    }

    out->fmt = fmt;
    out->data = data;
    out->data_bytes = data_bytes;
    return 0;
}

static int save_wav_pcm(const char *path, const WavFile *w) {
    FILE *f = fopen(path, "wb");
    if (!f) { perror("fopen"); return -1; }

    // We will write a minimal WAV: RIFF + fmt (16 bytes) + data
    const uint32_t fmt_chunk_size = 16;
    const uint32_t data_chunk_size = w->data_bytes;

    // RIFF size = 4 ("WAVE") + (8+fmt) + (8+data)
    uint32_t riff_size = 4 + (8 + fmt_chunk_size) + (8 + data_chunk_size);

    RiffHeader rh;
    memcpy(rh.riff_id, "RIFF", 4);
    rh.riff_size = riff_size;
    memcpy(rh.wave_id, "WAVE", 4);

    if (write_exact(f, &rh, sizeof(rh)) != 0) { fclose(f); return -1; }

    ChunkHeader fmt_h;
    memcpy(fmt_h.id, "fmt ", 4);
    fmt_h.size = fmt_chunk_size;

    if (write_exact(f, &fmt_h, sizeof(fmt_h)) != 0) { fclose(f); return -1; }
    if (write_exact(f, &w->fmt, fmt_chunk_size) != 0) { fclose(f); return -1; }

    ChunkHeader data_h;
    memcpy(data_h.id, "data", 4);
    data_h.size = data_chunk_size;

    if (write_exact(f, &data_h, sizeof(data_h)) != 0) { fclose(f); return -1; }
    if (write_exact(f, w->data, data_chunk_size) != 0) { fclose(f); return -1; }

    // pad if needed (minimal writers often skip; but we’ll be correct)
    if (data_chunk_size & 1) {
        uint8_t pad = 0;
        if (write_exact(f, &pad, 1) != 0) { fclose(f); return -1; }
    }

    fclose(f);
    return 0;
}

static void reverse_frames_pcm16(WavFile *w) {
    if (w->fmt.bits_per_sample != 16) {
        fprintf(stderr, "reverse_frames_pcm16: requires 16-bit PCM\n");
        return;
    }

    const uint32_t frame_bytes = w->fmt.block_align;  // channels * 2
    if (frame_bytes == 0 || (w->data_bytes % frame_bytes) != 0) {
        fprintf(stderr, "reverse_frames_pcm16: bad block_align/data_bytes\n");
        return;
    }

    uint32_t frames = w->data_bytes / frame_bytes;
    uint8_t *tmp = (uint8_t*)malloc(frame_bytes);
    if (!tmp) return;

    for (uint32_t i = 0; i < frames / 2; i++) {
        uint8_t *a = w->data + i * frame_bytes;
        uint8_t *b = w->data + (frames - 1 - i) * frame_bytes;
        memcpy(tmp, a, frame_bytes);
        memcpy(a, b, frame_bytes);
        memcpy(b, tmp, frame_bytes);
    }

    free(tmp);
}

static inline float princargf(float x) {
    // wrap to (-pi, pi]
    x = fmodf(x + (float)M_PI, 2.0f * (float)M_PI);
    if (x < 0) x += 2.0f * (float)M_PI;
    return x - (float)M_PI;
}

static inline float fast_atan2f(float y, float x) {
    return atan2f(y, x);
}

static void phase_vocoder(WavFile *w) {
    if (w->fmt.bits_per_sample != 16 || w->fmt.audio_format != 1) {
        fprintf(stderr, "phase_vocoder: requires 16-bit PCM WAV\n");
        return;
    }

    const int C = (int)w->fmt.num_channels;
    const int sr = (int)w->fmt.sample_rate;
    const uint32_t frame_bytes = w->fmt.block_align; // C * 2
    if (C <= 0 || frame_bytes != (uint32_t)(C * 2)) {
        fprintf(stderr, "phase_vocoder: unexpected block_align/channels\n");
        return;
    }
    if ((w->data_bytes % frame_bytes) != 0) {
        fprintf(stderr, "phase_vocoder: data_bytes not aligned to frames\n");
        return;
    }

    const int L = (int)(w->data_bytes / frame_bytes); // samples per channel
    const int16_t *pcm = (const int16_t*)w->data;

    // --------------------------
    // Phase vocoder parameters
    // --------------------------
    const int N  = 1024;   // FFT / window size
    const int Ha = 256;    // analysis hop (75% overlap)
    const float time_stretch = 0.50f; // >1.0 = longer, <1.0 = shorter
    const int Hs = (int)lroundf((float)Ha * time_stretch); // synthesis hop

    if (N <= 0 || Ha <= 0 || Hs <= 0 || Ha > N) {
        fprintf(stderr, "phase_vocoder: bad N/Ha/Hs\n");
        return;
    }

    const int K = N/2 + 1;
    const int num_frames = stft_num_frames(L, N, Ha);
    const int out_L = (num_frames - 1) * Hs + N; // output samples per channel

    // window (periodic Hann like your STFT)
    float *win = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    if (!win) return;
    hann_window(win, N);

    // FFTW buffers/plans
    float *time_buf = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    float *ifft_buf = (float*)fftwf_malloc(sizeof(float) * (size_t)N);
    fftwf_complex *X = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)K);
    fftwf_complex *Y = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)K);
    if (!time_buf || !ifft_buf || !X || !Y) {
        fftwf_free(win);
        if (time_buf) fftwf_free(time_buf);
        if (ifft_buf) fftwf_free(ifft_buf);
        if (X) fftwf_free(X);
        if (Y) fftwf_free(Y);
        return;
    }

    fftwf_plan p_r2c = fftwf_plan_dft_r2c_1d(N, time_buf, X, FFTW_MEASURE);
    fftwf_plan p_c2r = fftwf_plan_dft_c2r_1d(N, Y, ifft_buf, FFTW_MEASURE);
    if (!p_r2c || !p_c2r) {
        if (p_r2c) fftwf_destroy_plan(p_r2c);
        if (p_c2r) fftwf_destroy_plan(p_c2r);
        fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
        return;
    }

    // Output buffers per channel
    float *out = (float*)calloc((size_t)out_L * (size_t)C, sizeof(float));
    float *norm = (float*)calloc((size_t)out_L, sizeof(float)); // same for all channels
    if (!out || !norm) {
        free(out); free(norm);
        fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
        fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
        return;
    }

    // Phase vocoder state per channel per bin
    float *prev_phase = (float*)calloc((size_t)C * (size_t)K, sizeof(float));
    float *sum_phase  = (float*)calloc((size_t)C * (size_t)K, sizeof(float));
    if (!prev_phase || !sum_phase) {
        free(prev_phase); free(sum_phase);
        free(out); free(norm);
        fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
        fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
        return;
    }

    // Precompute expected phase advance per bin for analysis hop
    float *omega = (float*)malloc(sizeof(float) * (size_t)K);
    if (!omega) {
        free(prev_phase); free(sum_phase);
        free(out); free(norm);
        fftwf_destroy_plan(p_r2c); fftwf_destroy_plan(p_c2r);
        fftwf_free(win); fftwf_free(time_buf); fftwf_free(ifft_buf); fftwf_free(X); fftwf_free(Y);
        return;
    }
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
        for (int f = 0; f < num_frames; f++) {
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
        return;
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
}

static void print_info(const WavFile *w) {
    printf("WAV info:\n");
    printf("  format: PCM (1)\n");
    printf("  channels: %u\n", w->fmt.num_channels);
    printf("  sample_rate: %u\n", w->fmt.sample_rate);
    printf("  bits_per_sample: %u\n", w->fmt.bits_per_sample);
    printf("  block_align: %u\n", w->fmt.block_align);
    printf("  byte_rate: %u\n", w->fmt.byte_rate);
    printf("  data_bytes: %u\n", w->data_bytes);
    printf("  approx duration: %.3f sec\n",
           (w->fmt.byte_rate ? (double)w->data_bytes / (double)w->fmt.byte_rate : 0.0));
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s in.wav out.wav\n", argv[0]);
        return 1;
    }

    WavFile w;
    if (load_wav_pcm(argv[1], &w) != 0) {
        fprintf(stderr, "Failed to load wav\n");
        return 1;
    }

    print_info(&w);

    // --- EDIT HERE ---
    phase_vocoder(&w);

    if (save_wav_pcm(argv[2], &w) != 0) {
        fprintf(stderr, "Failed to save wav\n");
        free_wav(&w);
        return 1;
    }
    /* Resource handling */
    free_wav(&w);
    printf("Wrote edited file: %s\n", argv[2]);
    return 0;
}
