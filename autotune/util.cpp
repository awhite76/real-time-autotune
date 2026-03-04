#include "util.hpp"

#include <fstream>
#include <stdexcept>
#include <cstring>
#include <vector>

StereoWavI16 loadStereoWav_i16(const std::string &filename)
{
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file)
        throw std::runtime_error("Failed to open file");

    // --- RIFF ---
    char riff[4];
    file.read(riff, 4);
    if (std::strncmp(riff, "RIFF", 4) != 0)
        throw std::runtime_error("Not a RIFF file");

    file.ignore(4);

    char wave[4];
    file.read(wave, 4);
    if (std::strncmp(wave, "WAVE", 4) != 0)
        throw std::runtime_error("Not a WAVE file");

    uint16_t audioFormat = 0;
    uint16_t numChannels = 0;
    uint32_t sampleRate = 0;
    uint16_t bitsPerSample = 0;
    uint32_t dataSize = 0;
    std::streampos dataPos = 0;

    while (file)
    {
        char chunkId[4];
        uint32_t chunkSize = 0;

        file.read(chunkId, 4);
        file.read(reinterpret_cast<char *>(&chunkSize), 4);
        if (!file)
            break;

        if (std::strncmp(chunkId, "fmt ", 4) == 0)
        {
            if (chunkSize < 16)
                throw std::runtime_error("Invalid fmt chunk");

            file.read(reinterpret_cast<char *>(&audioFormat), 2);
            file.read(reinterpret_cast<char *>(&numChannels), 2);
            file.read(reinterpret_cast<char *>(&sampleRate), 4);

            file.ignore(6);
            file.read(reinterpret_cast<char *>(&bitsPerSample), 2);

            // Skip any extra fmt bytes
            if (chunkSize > 16)
                file.ignore(static_cast<std::streamsize>(chunkSize - 16));
        }
        else if (std::strncmp(chunkId, "data", 4) == 0)
        {
            dataSize = chunkSize;
            dataPos = file.tellg();
            break;
        }
        else
        {
            file.ignore(static_cast<std::streamsize>(chunkSize));
        }
    }

    if (audioFormat != 1)
        throw std::runtime_error("Only PCM supported");
    if (numChannels != 2)
        throw std::runtime_error("Not stereo");
    if (bitsPerSample != 16)
        throw std::runtime_error("Only 16-bit supported");
    if (dataSize == 0)
        throw std::runtime_error("Missing/empty data chunk");
    if (dataSize % 4 != 0)
        throw std::runtime_error("Invalid stereo data size");

    const size_t frames = dataSize / 4;

    file.seekg(dataPos);
    if (!file)
        throw std::runtime_error("Failed to seek to data");

    std::vector<int16_t> interleaved(frames * 2);
    file.read(reinterpret_cast<char *>(interleaved.data()),
              static_cast<std::streamsize>(dataSize));
    if (!file)
        throw std::runtime_error("Failed to read sample data");

    StereoWavI16 out;
    out.sampleRate = sampleRate;
    out.left.resize(frames);
    out.right.resize(frames);

    for (size_t i = 0; i < frames; ++i)
    {
        out.left[i] = interleaved[2 * i];
        out.right[i] = interleaved[2 * i + 1];
    }

    return out;
}

void deinterleave_stereo_i16(const int16_t *interleavedLR,
                             int16_t *left,
                             int16_t *right,
                             int frames)
{
    for (int i = 0; i < frames; ++i)
    {
        left[i] = interleavedLR[2 * i + 0];
        right[i] = interleavedLR[2 * i + 1];
    }
}

// Minimal WAV writer: 16-bit PCM stereo interleaved
static void writeStereoWav_i16_interleaved(
    const std::string &filename,
    uint32_t sampleRate,
    const int16_t *interleavedLR,
    uint32_t frames)
{
    if (!interleavedLR)
        throw std::runtime_error("writeStereoWav: null buffer");

    const uint16_t audioFormat = 1; // PCM
    const uint16_t numChannels = 2;
    const uint16_t bitsPerSample = 16;

    const uint32_t blockAlign = numChannels * (bitsPerSample / 8); // 4
    const uint32_t byteRate = sampleRate * blockAlign;
    const uint32_t dataSize = frames * blockAlign;

    const uint32_t riffSize = 4 /*WAVE*/ +
                              8 + 16 /*fmt*/ +
                              8 + dataSize /*data*/;

    std::ofstream out(filename, std::ios::binary);
    if (!out)
        throw std::runtime_error("Failed to open output wav: " + filename);

    // RIFF header
    out.write("RIFF", 4);
    out.write(reinterpret_cast<const char *>(&riffSize), 4);
    out.write("WAVE", 4);

    // fmt chunk
    out.write("fmt ", 4);
    uint32_t fmtSize = 16;
    out.write(reinterpret_cast<const char *>(&fmtSize), 4);
    out.write(reinterpret_cast<const char *>(&audioFormat), 2);
    out.write(reinterpret_cast<const char *>(&numChannels), 2);
    out.write(reinterpret_cast<const char *>(&sampleRate), 4);
    out.write(reinterpret_cast<const char *>(&byteRate), 4);
    out.write(reinterpret_cast<const char *>(&blockAlign), 2);
    out.write(reinterpret_cast<const char *>(&bitsPerSample), 2);

    // data chunk
    out.write("data", 4);
    out.write(reinterpret_cast<const char *>(&dataSize), 4);
    out.write(reinterpret_cast<const char *>(interleavedLR), dataSize);

    if (!out)
        throw std::runtime_error("Failed while writing wav: " + filename);
}

static std::vector<int16_t> interleaveStereo(const StereoWavI16 &wav)
{
    if (wav.left.size() != wav.right.size())
        throw std::runtime_error("interleaveStereo: L/R size mismatch");
    std::vector<int16_t> out(wav.left.size() * 2);
    for (size_t i = 0; i < wav.left.size(); ++i)
    {
        out[2 * i + 0] = wav.left[i];
        out[2 * i + 1] = wav.right[i];
    }
    return out;
}