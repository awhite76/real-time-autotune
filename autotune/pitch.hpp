#ifndef pitch_h
#define pitch_h

#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>

using namespace std;

struct PitchSeries
{
    uint32_t sampleRate = 0;
    uint32_t windowFrames = 0; // frames per chunk (per channel)
    float windowMs = 0.0f;

    // pitch in Hz per chunk (same indexing for left/right)
    vector<float> leftHz;
    vector<float> rightHz;

    // optional: confidence per chunk (if you want it)
    vector<float> leftConf;
    vector<float> rightConf;

    // quick helpers
    inline size_t size() const { return leftHz.size(); }

    // map time (ms) -> chunk index
    inline size_t indexForMs(float tMs) const
    {
        if (windowMs <= 0.0f)
            return 0;
        size_t idx = static_cast<size_t>(tMs / windowMs);
        return (idx >= size()) ? (size() ? size() - 1 : 0) : idx;
    }
};

class Yin
{

public:
    Yin();
    Yin(int bufferSize);
    void initialize(int bufferSize);
    float getPitch(int16_t *buffer);
    float getProbability();

private:
    float parabolicInterpolation(int tauEstimate);
    int absoluteThreshold();
    void cumulativeMeanNormalizedDifference();
    void difference(int16_t *buffer);

    double threshold;
    int bufferSize;
    int halfBufferSize;
    float sampleRate;
    float *yinBuffer;
    float probability;
};

#endif