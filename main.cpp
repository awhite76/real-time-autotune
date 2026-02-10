#include <iostream>
#include <alsa/asoundlib.h>
#include <cstdint>
#include <cstring>
#include <string>

#define T_MS 15
#define SAMPLE_RATE 48000
#define CHANNELS 1
#define PERIOD_FRAMES (SAMPLE_RATE * T_MS / 1000)
#define BUFFER_FRAMES (PERIOD_FRAMES * 2)

static bool set_hw_params(snd_pcm_t *handle, snd_pcm_stream_t stream)
{
    snd_pcm_hw_params_t *hw = nullptr;
    snd_pcm_hw_params_alloca(&hw);

    int err = 0;
    if ((err = snd_pcm_hw_params_any(handle, hw)) < 0)
    {
        cerr << "hw_params_any: " << snd_strerror(err) << "\n";
        return false;
    }

    // Interleaved, S16_LE, single channel, 48k
    if ((err = snd_pcm_hw_params_set_access(handle, hw, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0)
    {
        cerr << "set_access: " << snd_strerror(err) << "\n";
        return false;
    }
    if ((err = snd_pcm_hw_params_set_format(handle, hw, SND_PCM_FORMAT_S16_LE)) < 0)
    {
        cerr << "set_format: " << snd_strerror(err) << "\n";
        return false;
    }
    if ((err = snd_pcm_hw_params_set_channels(handle, hw, CHANNELS)) < 0)
    {
        cerr << "set_channels: " << snd_strerror(err) << "\n";
        return false;
    }

    unsigned int rate = SAMPLE_RATE;
    if ((err = snd_pcm_hw_params_set_rate_near(handle, hw, &rate, nullptr)) < 0)
    {
        cerr << "set_rate_near: " << snd_strerror(err) << "\n";
        return false;
    }
    if (rate != SAMPLE_RATE)
    {
        cerr << "Warning: device set rate to " << rate << " (requested " << SAMPLE_RATE << ")\n";
    }

    // Period/buffer sizes
    snd_pcm_uframes_t period = PERIOD_FRAMES;
    if ((err = snd_pcm_hw_params_set_period_size_near(handle, hw, &period, nullptr)) < 0)
    {
        cerr << "set_period_size_near: " << snd_strerror(err) << "\n";
        return false;
    }

    snd_pcm_uframes_t buffer = BUFFER_FRAMES;
    if ((err = snd_pcm_hw_params_set_buffer_size_near(handle, hw, &buffer)) < 0)
    {
        cerr << "set_buffer_size_near: " << snd_strerror(err) << "\n";
        return false;
    }

    if ((err = snd_pcm_hw_params(handle, hw)) < 0)
    {
        cerr << "hw_params: " << snd_strerror(err) << "\n";
        return false;
    }

    // Software params: start ASAP, low latency behavior
    snd_pcm_sw_params_t *sw = nullptr;
    snd_pcm_sw_params_alloca(&sw);

    if ((err = snd_pcm_sw_params_current(handle, sw)) < 0)
    {
        cerr << "sw_params_current: " << snd_strerror(err) << "\n";
        return false;
    }

    // Start playback when at least one period is available
    if (stream == SND_PCM_STREAM_PLAYBACK)
    {
        if ((err = snd_pcm_sw_params_set_start_threshold(handle, sw, PERIOD_FRAMES)) < 0)
        {
            cerr << "set_start_threshold: " << snd_strerror(err) << "\n";
            return false;
        }
    }

    // Wake up when at least a period can be processed
    if ((err = snd_pcm_sw_params_set_avail_min(handle, sw, PERIOD_FRAMES)) < 0)
    {
        cerr << "set_avail_min: " << snd_strerror(err) << "\n";
        return false;
    }

    if ((err = snd_pcm_sw_params(handle, sw)) < 0)
    {
        cerr << "sw_params: " << snd_strerror(err) << "\n";
        return false;
    }

    return true;
}

int main(int argc, char **argv)
{
    // Usage: ./main <capture_dev> <playback_dev>
    string cap_dev = (argc > 1) ? argv[1] : "hw:0,0";
    string pb_dev = (argc > 2) ? argv[2] : "hw:0,0";

    snd_pcm_t *capture_handle = nullptr;
    snd_pcm_t *playback_handle = nullptr;

    int err = 0;

    // Open devices
    if ((err = snd_pcm_open(&capture_handle, cap_dev.c_str(), SND_PCM_STREAM_CAPTURE, 0)) < 0)
    {
        cerr << "snd_pcm_open CAPTURE (" << cap_dev << "): " << snd_strerror(err) << "\n";
        return 1;
    }
    if ((err = snd_pcm_open(&playback_handle, pb_dev.c_str(), SND_PCM_STREAM_PLAYBACK, 0)) < 0)
    {
        cerr << "snd_pcm_open PLAYBACK (" << pb_dev << "): " << snd_strerror(err) << "\n";
        snd_pcm_close(capture_handle);
        return 1;
    }

    // Configure HW/SW params
    if (!set_hw_params(capture_handle, SND_PCM_STREAM_CAPTURE) ||
        !set_hw_params(playback_handle, SND_PCM_STREAM_PLAYBACK))
    {
        snd_pcm_close(playback_handle);
        snd_pcm_close(capture_handle);
        return 1;
    }

    // Prepare devices
    if ((err = snd_pcm_prepare(capture_handle)) < 0)
    {
        cerr << "prepare capture: " << snd_strerror(err) << "\n";
        snd_pcm_close(playback_handle);
        snd_pcm_close(capture_handle);
        return 1;
    }
    if ((err = snd_pcm_prepare(playback_handle)) < 0)
    {
        cerr << "prepare playback: " << snd_strerror(err) << "\n";
        snd_pcm_close(playback_handle);
        snd_pcm_close(capture_handle);
        return 1;
    }

    int16_t buffer[PERIOD_FRAMES * CHANNELS];
    memset(buffer, 0, sizeof(buffer));
    snd_pcm_writei(playback_handle, buffer, PERIOD_FRAMES);

    cerr << "Looping " << T_MS << "ms chunks: "
         << PERIOD_FRAMES << " frames @ " << SAMPLE_RATE << " Hz, "
         << CHANNELS << " ch\n"
         << "Capture: " << cap_dev << "  Playback: " << pb_dev << "\n";

    // Main loop: read 15ms, write 15ms
    while (true)
    {
        // Read exactly PERIOD_FRAMES frames
        snd_pcm_sframes_t rcvd = 0;
        while (rcvd < PERIOD_FRAMES)
        {
            snd_pcm_sframes_t r = snd_pcm_readi(
                capture_handle,
                buffer + rcvd * CHANNELS,
                PERIOD_FRAMES - rcvd);
            if (r < 0)
            {
                cerr << "capture read failed: " << snd_strerror((int)r) << "\n";
                goto out;
            }
            rcvd += r;
        }

        // TODO: process buffer here

        // Write exactly PERIOD_FRAMES frames
        snd_pcm_sframes_t sent = 0;
        while (sent < PERIOD_FRAMES)
        {
            snd_pcm_sframes_t w = snd_pcm_writei(
                playback_handle,
                buffer + sent * CHANNELS,
                PERIOD_FRAMES - sent);
            if (w < 0)
            {
                cerr << "playback write failed: " << snd_strerror((int)w) << "\n";
                goto out;
            }
            sent += w;
        }
    }

out:
    snd_pcm_drain(playback_handle);
    snd_pcm_close(playback_handle);
    snd_pcm_close(capture_handle);
    return 0;
}
