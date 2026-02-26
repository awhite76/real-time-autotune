#include <iostream>
#include <alsa/asoundlib.h>
#include <cstdint>
#include <cstring>
#include <string>
#include <algorithm>

#define T_MS 10
#define SAMPLE_RATE 48000
#define CHANNELS 2
#define PERIOD_FRAMES (SAMPLE_RATE * T_MS / 1000)
#define BUFFER_FRAMES (PERIOD_FRAMES * CHANNELS * 2)

using namespace std;

bool set_hw_params(snd_pcm_t *handle, snd_pcm_stream_t stream);
int xrun_recover(snd_pcm_t *handle, int err);
