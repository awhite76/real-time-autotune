#ifndef TS_H
#define TS_H

#include <speex/speex_resampler.h>
#include <cstdint>

struct TimeStretchResampler
{
    SpeexResamplerState *st = nullptr;
    unsigned int sampleRate = 0;
    float last_s = -1.0f; // cache last ratio to avoid reconfig every call
};

/**
 * Initialize SpeexDSP resampler for interleaved multi-channel audio.
 *
 * The number of channels is taken from CHANNELS (main.hpp).
 *
 * @param r          Resampler state
 * @param sampleRate Base sample rate (Hz)
 * @param quality    0..10 (5 is a good realtime compromise)
 * @return true on success
 */
bool time_stretch_init(TimeStretchResampler &r,
                       unsigned int sampleRate,
                       int quality = 5);

/**
 * Destroy resampler state (safe to call multiple times).
 */
void time_stretch_destroy(TimeStretchResampler &r);

/**
 * Resample interleaved int16 audio by factor s.
 *
 * NOTE: This is resampling (changes pitch + duration). In a PV pipeline, you
 * typically time-scale by 1/s then resample by s to achieve pitch shift with
 * near-constant duration.
 *
 * Buffers are interleaved with CHANNELS channels:
 *   [ch0, ch1, ..., chN-1, ch0, ch1, ..., chN-1, ...]
 *
 * @param r           Resampler state
 * @param input       Interleaved int16 input
 * @param inFrames    Input frames per channel
 * @param output      Interleaved int16 output buffer
 * @param outCapacity Output capacity in frames per channel
 * @param s           Resample ratio (>0)
 * @return output frames written per channel (0 on error)
 */
int time_stretch_process(TimeStretchResampler &r,
                         const int16_t *input,
                         int inFrames,
                         int16_t *output,
                         int outCapacity,
                         float s);

#endif