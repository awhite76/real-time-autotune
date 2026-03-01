#include "time_stretch.hpp"
#include <cmath>
#include "main.hpp"

bool time_stretch_init(TimeStretchResampler &r,
                       unsigned int sampleRate,
                       int quality)
{
    int err = 0;

    r.sampleRate = sampleRate;
    r.last_s = -1.0f;

    r.st = speex_resampler_init(
        CHANNELS,   // channels
        sampleRate, // in rate
        sampleRate, // out rate (we'll change ratio dynamically)
        quality,    // 0..10
        &err);

    return (err == RESAMPLER_ERR_SUCCESS);
}

void time_stretch_destroy(TimeStretchResampler &r)
{
    if (r.st)
        speex_resampler_destroy(r.st);

    r.st = nullptr;
    r.sampleRate = 0;
    r.last_s = -1.0f;
}

int time_stretch_process(TimeStretchResampler &r,
                         const int16_t *input,
                         int inFrames, // frames per channel
                         int16_t *output,
                         int outCapacity, // frames per channel
                         float s)
{
    if (!r.st || !input || !output ||
        inFrames <= 0 || outCapacity <= 0 || s <= 0.0f)
        return 0;


    cerr << "Out L is " << inFrames << "\n";
    cerr << "s is " << s << "\n";
    cerr << "frames per channel is " << outCapacity << "\n";
    // Update ratio only if changed
    if (s != r.last_s)
    {

        cerr << "Before set rate frac" << "\n";
        const int den = 100;
        int num = (int)std::lround((double)s * (double)den);
        if (num < 1)
            num = 1;

        speex_resampler_set_rate_frac(
            r.st,
            num,
            den,
            r.sampleRate,
            r.sampleRate);

        r.last_s = s;
    }

    // IMPORTANT:
    // inFrames/outCapacity are FRAMES PER CHANNEL
    // Buffers are interleaved [L R L R ...]
    spx_uint32_t inLen = (spx_uint32_t)inFrames;
    spx_uint32_t outLen = (spx_uint32_t)outCapacity;

    cerr << "Before resampler process interleavened......" << "\n";
    cerr << "type casted inLen" << inLen << "\n";
    cerr << "type casted outLen" << outLen << "\n";
    int err = speex_resampler_process_interleaved_int(
        r.st,
        input,
        &inLen,
        output,
        &outLen);

    if (err != RESAMPLER_ERR_SUCCESS)
        return 0;

    return (int)outLen; // frames per channel
}
