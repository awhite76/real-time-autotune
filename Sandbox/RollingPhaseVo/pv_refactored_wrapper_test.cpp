#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

//g++ -O2 -std=c++17 pv_refactored_wrapper_test.cpp pv_refactor.cpp -lfftw3f -lm -o pv_test
//   run : ./pv_test input.wav output.wav 2.0 256
// Your PV API
#include "pv_refactor.hpp"

// ---------------- WAV (PCM16 mono) minimal reader/writer ----------------

#pragma pack(push, 1)
struct WavHeader {
    char     riff[4];        // "RIFF"
    uint32_t riffSize;       // fileSize - 8
    char     wave[4];        // "WAVE"

    char     fmt[4];         // "fmt "
    uint32_t fmtSize;        // 16 for PCM
    uint16_t audioFormat;    // 1 = PCM
    uint16_t numChannels;    // 1
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;  // 16

    char     data[4];        // "data"
    uint32_t dataSize;       // numSamples * numChannels * bits/8
};
#pragma pack(pop)

static uint32_t read_u32_le(FILE* f) {
    uint8_t b[4];
    if (fread(b, 1, 4, f) != 4) return 0;
    return (uint32_t)b[0] | ((uint32_t)b[1] << 8) | ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static uint16_t read_u16_le(FILE* f) {
    uint8_t b[2];
    if (fread(b, 1, 2, f) != 2) return 0;
    return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
}

static bool read_wav_pcm16_mono(const char* path, std::vector<int16_t>& samples, uint32_t& sampleRate) {
    samples.clear();
    sampleRate = 0;

    FILE* f = std::fopen(path, "rb");
    if (!f) return false;

    char riff[4], wave[4];
    if (fread(riff, 1, 4, f) != 4) { fclose(f); return false; }
    (void)read_u32_le(f); // riffSize
    if (fread(wave, 1, 4, f) != 4) { fclose(f); return false; }

    if (std::memcmp(riff, "RIFF", 4) != 0 || std::memcmp(wave, "WAVE", 4) != 0) {
        fclose(f);
        return false;
    }

    bool got_fmt = false;
    bool got_data = false;

    // fmt fields we care about
    uint16_t audioFormat = 0;
    uint16_t numChannels = 0;
    uint16_t bitsPerSample = 0;

    std::vector<uint8_t> data_bytes;

    // Scan chunks
    while (!got_data) {
        char chunkId[4];
        if (fread(chunkId, 1, 4, f) != 4) break;
        uint32_t chunkSize = read_u32_le(f);

        // Chunks are word-aligned: if chunkSize is odd, there is 1 pad byte.
        const uint32_t paddedSize = (chunkSize + 1u) & ~1u;

        if (std::memcmp(chunkId, "fmt ", 4) == 0) {
            // PCM fmt chunk is at least 16 bytes
            audioFormat   = read_u16_le(f);
            numChannels   = read_u16_le(f);
            sampleRate    = read_u32_le(f);
            (void)read_u32_le(f); // byteRate
            (void)read_u16_le(f); // blockAlign
            bitsPerSample = read_u16_le(f);

            // Skip any remaining fmt bytes
            const uint32_t consumed = 16;
            if (chunkSize > consumed) {
                fseek(f, (long)(chunkSize - consumed), SEEK_CUR);
            }

            got_fmt = true;
        }
        else if (std::memcmp(chunkId, "data", 4) == 0) {
            data_bytes.resize(chunkSize);
            if (chunkSize > 0) {
                if (fread(data_bytes.data(), 1, chunkSize, f) != chunkSize) {
                    fclose(f);
                    return false;
                }
            }
            got_data = true;
        }
        else {
            // Skip unhandled chunk
            fseek(f, (long)paddedSize, SEEK_CUR);
        }

        // If we processed exactly chunkSize, we still need to skip pad byte if odd
        // (Handled by paddedSize when skipping, but fmt/data branches consumed directly.)
        if ((std::memcmp(chunkId, "fmt ", 4) == 0) || (std::memcmp(chunkId, "data", 4) == 0)) {
            if (chunkSize & 1u) fseek(f, 1, SEEK_CUR);
        }
    }

    fclose(f);

    if (!got_fmt || !got_data) return false;

    // Only accept PCM format=1 for this minimal loader
    if (audioFormat != 1) return false;
    if (numChannels != 1) return false;
    if (bitsPerSample != 16) return false;

    // Convert bytes to int16_t samples
    if (data_bytes.size() % 2 != 0) return false;
    const size_t nSamp = data_bytes.size() / 2;
    samples.resize(nSamp);
    for (size_t i = 0; i < nSamp; i++) {
        uint16_t lo = data_bytes[2*i + 0];
        uint16_t hi = data_bytes[2*i + 1];
        samples[i] = (int16_t)(lo | (hi << 8)); // little-endian
    }

    return true;
}

static bool write_wav_pcm16_mono(const char* path, const std::vector<int16_t>& samples, uint32_t sampleRate) {
    FILE* f = std::fopen(path, "wb");
    if (!f) return false;

    WavHeader h{};
    std::memcpy(h.riff, "RIFF", 4);
    std::memcpy(h.wave, "WAVE", 4);
    std::memcpy(h.fmt,  "fmt ", 4);
    std::memcpy(h.data, "data", 4);

    h.fmtSize       = 16;
    h.audioFormat   = 1;
    h.numChannels   = 1;
    h.sampleRate    = sampleRate;
    h.bitsPerSample = 16;
    h.blockAlign    = (h.numChannels * h.bitsPerSample) / 8;
    h.byteRate      = h.sampleRate * h.blockAlign;

    h.dataSize      = (uint32_t)(samples.size() * sizeof(int16_t));
    h.riffSize      = 4 + (8 + h.fmtSize) + (8 + h.dataSize);

    if (std::fwrite(&h, sizeof(h), 1, f) != 1) { std::fclose(f); return false; }
    if (!samples.empty()) {
        if (std::fwrite(samples.data(), sizeof(int16_t), samples.size(), f) != samples.size()) {
            std::fclose(f);
            return false;
        }
    }

    std::fclose(f);
    return true;
}

// ---------------- Main: stream simulation ----------------

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr,
            "Usage: %s <in_mono16.wav> <out_mono16.wav> <time_stretch> [chunk]\n"
            "Example: %s in.wav out.wav 2.0 256\n",
            argv[0], argv[0]);
        return 1;
    }

    const char* inPath  = argv[1];
    const char* outPath = argv[2];
    const float stretch = std::strtof(argv[3], nullptr);
    const size_t chunk  = (argc >= 5) ? (size_t)std::strtoul(argv[4], nullptr, 10) : 256;

    if (!(stretch > 0.0f)) {
        std::fprintf(stderr, "time_stretch must be > 0\n");
        return 1;
    }

    // 1) Load WAV into memory
    std::vector<int16_t> in;
    uint32_t fs = 0;
    if (!read_wav_pcm16_mono(inPath, in, fs)) {
        std::fprintf(stderr, "Failed to read WAV (must be PCM16 mono): %s\n", inPath);
        return 1;
    }

    // 2) Setup vocoder
    PhaseVocoder pv = (PhaseVocoder)std::calloc(1, sizeof(*pv));
    if (!pv) return 1;

    if (setup_vocoder(pv) != 0) {
        std::fprintf(stderr, "setup_vocoder failed\n");
        std::free(pv);
        return 1;
    }

    pv->time_stretch = stretch;

    // 3) Stream simulation: push chunks, process ready frames, collect output
    std::vector<int16_t> out;
    out.reserve((size_t)std::ceil((double)in.size() * (double)stretch) + WINDOW_SIZE);

    std::vector<int16_t> tmp(8192);

    size_t pos = 0;
    while (pos < in.size()) {
        const size_t n = std::min(chunk, in.size() - pos);

        size_t wrote = pv_push_input(pv, in.data() + pos, n);
        pos += wrote;

        // Pull whatever is currently ready
        size_t produced = pv_process_ready(pv, tmp.data(), tmp.size());
        out.insert(out.end(), tmp.begin(), tmp.begin() + produced);

        // If the ring was full and we couldn't write all n, process and retry remainder
        while (wrote < n) {
            produced = pv_process_ready(pv, tmp.data(), tmp.size());
            out.insert(out.end(), tmp.begin(), tmp.begin() + produced);

            const size_t remaining = n - wrote;
            const size_t more = pv_push_input(pv, in.data() + pos, remaining);
            pos += more;
            wrote += more;

            if (more == 0) break;
        }
    }

    // 4) Flush tail: feed zeros so final frames can run
    std::vector<int16_t> zeros(chunk, 0);
    size_t flush = WINDOW_SIZE + ANALYSIS_HOP;
    while (flush > 0) {
        const size_t n = std::min(chunk, flush);
        pv_push_input(pv, zeros.data(), n);
        flush -= n;

        size_t produced = pv_process_ready(pv, tmp.data(), tmp.size());
        out.insert(out.end(), tmp.begin(), tmp.begin() + produced);
    }

    // 5) Optional: drain a bit of overlap-add tail
    {
        std::vector<int16_t> tail(WINDOW_SIZE);
        pv_consume_output(pv, tail.data(), tail.size());
        out.insert(out.end(), tail.begin(), tail.end());
    }

    // 6) Cleanup
    cleanup_vocoder(pv);
    std::free(pv);

    // 7) Write output WAV
    if (!write_wav_pcm16_mono(outPath, out, fs)) {
        std::fprintf(stderr, "Failed to write WAV: %s\n", outPath);
        return 1;
    }

    std::fprintf(stderr, "Done. in=%zu samples, out=%zu samples, stretch=%.3f\n",
                 in.size(), out.size(), stretch);
    return 0;
}