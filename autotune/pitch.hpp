#ifndef pitch_h
#define pitch_h

#include <cstdlib>
#include <cstring>
#include <cstdint>

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