#ifndef UTIL_H
#define UTIL_H

#include <cstdint>
#include <string>
#include <vector>

struct StereoWavI16
{
    uint32_t sampleRate = 0;
    std::vector<int16_t> left;
    std::vector<int16_t> right;
};

// Loads stereo 16-bit PCM WAV.
// Throws std::runtime_error on failure.
StereoWavI16 loadStereoWav_i16(const std::string &filename);

// Deinterleave interleaved stereo LRLR... into L and R buffers.
void deinterleave_stereo_i16(const int16_t *interleavedLR,
                             int16_t *left,
                             int16_t *right,
                             int frames);

#endif