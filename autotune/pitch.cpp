#include "pitch.hpp"
#include "main.hpp"

void Yin::initialize(int yinBufferSize)
{
    bufferSize = yinBufferSize;
    sampleRate = SAMPLE_RATE;
    halfBufferSize = bufferSize / 2;
    threshold = 0.15f;
    probability = 0.0f;

    yinBuffer = (float *)malloc(sizeof(float) * halfBufferSize);

    for (int i = 0; i < halfBufferSize; i++)
        yinBuffer[i] = 0.0f;
}

Yin::Yin()
{
}

Yin::Yin(int yinBufferSize)
{
    initialize(yinBufferSize);
}

float Yin::getProbability()
{
    return probability;
}

float Yin::getPitch(int16_t *buffer)
{
    int tauEstimate = -1;
    float pitchInHertz = -1.0f;

    difference(buffer);
    cumulativeMeanNormalizedDifference();
    tauEstimate = absoluteThreshold();

    if (tauEstimate != -1)
        pitchInHertz = sampleRate / parabolicInterpolation(tauEstimate);

    return pitchInHertz;
}

float Yin::parabolicInterpolation(int tauEstimate)
{
    float betterTau;
    int x0;
    int x2;

    if (tauEstimate < 1)
        x0 = tauEstimate;
    else
        x0 = tauEstimate - 1;

    if (tauEstimate + 1 < halfBufferSize)
        x2 = tauEstimate + 1;
    else
        x2 = tauEstimate;

    if (x0 == tauEstimate)
    {
        betterTau = (yinBuffer[tauEstimate] <= yinBuffer[x2])
                        ? tauEstimate
                        : x2;
    }
    else if (x2 == tauEstimate)
    {
        betterTau = (yinBuffer[tauEstimate] <= yinBuffer[x0])
                        ? tauEstimate
                        : x0;
    }
    else
    {
        float s0 = yinBuffer[x0];
        float s1 = yinBuffer[tauEstimate];
        float s2 = yinBuffer[x2];

        betterTau = tauEstimate +
                    (s2 - s0) /
                        (2.0f * (2.0f * s1 - s2 - s0));
    }

    return betterTau;
}

void Yin::cumulativeMeanNormalizedDifference()
{
    yinBuffer[0] = 1.0f;
    float runningSum = 0.0f;

    for (int tau = 1; tau < halfBufferSize; tau++)
    {
        runningSum += yinBuffer[tau];

        if (runningSum > 0.0f)
            yinBuffer[tau] *= (float)tau / runningSum;
        else
            yinBuffer[tau] = 1.0f; // silence protection
    }
}

void Yin::difference(int16_t *buffer)
{
    // Clear buffer every run
    for (int tau = 0; tau < halfBufferSize; tau++)
        yinBuffer[tau] = 0.0f;

    for (int tau = 0; tau < halfBufferSize; tau++)
    {
        for (int index = 0; index < halfBufferSize; index++)
        {
            // use 32-bit intermediate to avoid overflow
            int32_t delta =
                (int32_t)buffer[index] -
                (int32_t)buffer[index + tau];

            float fdelta = (float)delta;
            yinBuffer[tau] += fdelta * fdelta;
        }
    }
}

int Yin::absoluteThreshold()
{
    int tau;

    for (tau = 2; tau < halfBufferSize; tau++)
    {
        if (yinBuffer[tau] < threshold)
        {
            while (tau + 1 < halfBufferSize &&
                   yinBuffer[tau + 1] < yinBuffer[tau])
            {
                tau++;
            }

            probability = 1.0f - yinBuffer[tau];
            break;
        }
    }

    if (tau == halfBufferSize || yinBuffer[tau] >= threshold)
    {
        tau = -1;
        probability = 0.0f;
    }

    return tau;
}