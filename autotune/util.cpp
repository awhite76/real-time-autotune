#include "util.hpp"

#include <fstream>
#include <stdexcept>
#include <cstring>

using namespace std;

StereoWavI16 loadStereoWav_i16(const string &filename)
{
    ifstream file(filename, ios::binary);
    if (!file)
        throw runtime_error("Failed to open file");

    // --- RIFF ---
    char riff[4];
    file.read(riff, 4);
    if (strncmp(riff, "RIFF", 4) != 0)
        throw runtime_error("Not a RIFF file");

    file.ignore(4);

    char wave[4];
    file.read(wave, 4);
    if (strncmp(wave, "WAVE", 4) != 0)
        throw runtime_error("Not a WAVE file");

    uint16_t audioFormat = 0;
    uint16_t numChannels = 0;
    uint32_t sampleRate = 0;
    uint16_t bitsPerSample = 0;
    uint32_t dataSize = 0;
    streampos dataPos = 0;

    while (file)
    {
        char chunkId[4];
        uint32_t chunkSize;

        file.read(chunkId, 4);
        file.read(reinterpret_cast<char *>(&chunkSize), 4);
        if (!file)
            break;

        if (strncmp(chunkId, "fmt ", 4) == 0)
        {
            file.read(reinterpret_cast<char *>(&audioFormat), 2);
            file.read(reinterpret_cast<char *>(&numChannels), 2);
            file.read(reinterpret_cast<char *>(&sampleRate), 4);

            file.ignore(6);
            file.read(reinterpret_cast<char *>(&bitsPerSample), 2);

            file.ignore(chunkSize - 16);
        }
        else if (strncmp(chunkId, "data", 4) == 0)
        {
            dataSize = chunkSize;
            dataPos = file.tellg();
            break;
        }
        else
        {
            file.ignore(chunkSize);
        }
    }

    if (audioFormat != 1)
        throw runtime_error("Only PCM supported");
    if (numChannels != 2)
        throw runtime_error("Not stereo");
    if (bitsPerSample != 16)
        throw runtime_error("Only 16-bit supported");
    if (dataSize % 4 != 0)
        throw runtime_error("Invalid stereo data size");

    size_t frames = dataSize / 4;

    file.seekg(dataPos);

    vector<int16_t> interleaved(frames * 2);
    file.read(reinterpret_cast<char *>(interleaved.data()), dataSize);

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