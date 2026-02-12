/* Boilerplate code from ChatGPT */

// wav_edit.c
// Build: gcc -O2 -Wall wave_edit.c -o wave_edit
// Run:   ./wave_edit in.wav out.wav

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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
    reverse_frames_pcm16(&w);

    if (save_wav_pcm(argv[2], &w) != 0) {
        fprintf(stderr, "Failed to save wav\n");
        free_wav(&w);
        return 1;
    }

    free_wav(&w);
    printf("Wrote edited file: %s\n", argv[2]);
    return 0;
}
