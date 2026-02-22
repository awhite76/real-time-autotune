#ifndef TS_H
#define TS_H

#include <speex/speex_resampler.h>
#include <cstdint>

struct TimeStretchResampler
{
    SpeexResamplerState *st = nullptr;
    unsigned int sampleRate = 0;

    // Optional: cache last ratio so we only update when it changes
    float last_s = -1.0f;
};

/**
 * Initialize a mono SpeexDSP resampler for dynamic ratio changes.
 * @param r          Resampler state
 * @param sampleRate Input/output base sample rate (Hz)
 * @param quality    0..10 (5 is a good realtime compromise)
 * @return true on success
 */
bool time_stretch_init(TimeStretchResampler &r,
                       unsigned int sampleRate,
                       int quality = 5);

/**
 * Destroy the resampler state (safe to call multiple times).
 */
void time_stretch_destroy(TimeStretchResampler &r);

/**
 * Resample by factor s (NOTE: this is resampling, not true time-stretch).
 *
 * s > 1.0 => output shorter & higher pitch
 * s < 1.0 => output longer  & lower pitch
 *
 * @param r           Resampler state
 * @param input       Mono int16 input samples
 * @param inFrames    Number of input frames (samples)
 * @param output      Output buffer (mono int16)
 * @param outCapacity Output capacity in frames (samples)
 * @param s           Ratio factor (>0)
 * @return number of output frames written (0 on error)
 */
int time_stretch_process(TimeStretchResampler &r,
                         const int16_t *input,
                         int inFrames,
                         int16_t *output,
                         int outCapacity,
                         float s);

#endif