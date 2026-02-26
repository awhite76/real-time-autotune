#include <iostream>
#include <alsa/asoundlib.h>
#include <cstdint>
#include <cstring>
#include <string>

/* Sampling data */
#define T_MS 20
#define SAMPLE_RATE 48000
#define CHANNELS 2
#define PERIOD_FRAMES (SAMPLE_RATE * T_MS / 1000) // Period frames represents how many samples per sampling period of one channel

/* DMA engine buffer */
#define BUFFER_FRAMES (PERIOD_FRAMES * CHANNELS) // Buffer frames represents size of DMA ring buffer

using namespace std;

bool set_hw_params(snd_pcm_t *handle, snd_pcm_stream_t stream);
int xrun_recover(snd_pcm_t *handle, int err);
