#include "main.hpp"

static int xrun_recover(snd_pcm_t *handle, int err)
{
    if (err == -EPIPE)
    { // XRUN
        cerr << "XRUN (-EPIPE). Preparing device...\n";
        err = snd_pcm_prepare(handle);
        return err;
    }
    if (err == -ESTRPIPE)
    { // suspend
        cerr << "Stream suspended (-ESTRPIPE). Resuming...\n";
        while ((err = snd_pcm_resume(handle)) == -EAGAIN)
        {
            snd_pcm_wait(handle, 100);
        }
        if (err < 0)
            err = snd_pcm_prepare(handle);
        return err;
    }
    return err;
}

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

    // Software params
    snd_pcm_sw_params_t *sw = nullptr;
    snd_pcm_sw_params_alloca(&sw);

    if ((err = snd_pcm_sw_params_current(handle, sw)) < 0)
    {
        cerr << "sw_params_current: " << snd_strerror(err) << "\n";
        return false;
    }

    // Start playback only after we have at least 2 periods queued
    if (stream == SND_PCM_STREAM_PLAYBACK)
    {
        if ((err = snd_pcm_sw_params_set_start_threshold(handle, sw, PERIOD_FRAMES * 2)) < 0)
        {
            cerr << "set_start_threshold: " << snd_strerror(err) << "\n";
            return false;
        }
    }

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

    /* Speaker settup */
    string cap_dev = (argc > 1) ? argv[1] : "hw:2,0";
    string pb_dev = (argc > 2) ? argv[2] : "hw:2,0";

    snd_pcm_t *capture_handle = nullptr;
    snd_pcm_t *playback_handle = nullptr;

    int err = 0;

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

    if (!set_hw_params(capture_handle, SND_PCM_STREAM_CAPTURE) ||
        !set_hw_params(playback_handle, SND_PCM_STREAM_PLAYBACK))
    {
        snd_pcm_close(playback_handle);
        snd_pcm_close(capture_handle);
        return 1;
    }

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

    // Prefill playback with 2 periods of silence so it won't underrun while capture blocks
    for (int i = 0; i < 2; i++)
    {
        snd_pcm_sframes_t w = snd_pcm_writei(playback_handle, buffer, PERIOD_FRAMES);
        if (w < 0)
        {
            w = xrun_recover(playback_handle, (int)w);
            if (w < 0)
            {
                cerr << "initial playback prefill failed: " << snd_strerror((int)w) << "\n";
                snd_pcm_drop(playback_handle);
                snd_pcm_close(playback_handle);
                snd_pcm_close(capture_handle);
                return 0;
            }
            i--; // retry this period after recover
        }
    }

    cerr << "Looping " << T_MS << "ms chunks: "
         << PERIOD_FRAMES << " frames @ " << SAMPLE_RATE << " Hz, "
         << CHANNELS << " ch\n"
         << "Capture: " << cap_dev << "  Playback: " << pb_dev << "\n";

    /* Phase Vo settup */

    float *time_buf; float *win; float *ifft_buf; float *omega; float *out;
    float *norm; int16_t *new_data; float *prev_phase; float *sum_phase;
    fftwf_complex *X; fftwf_complex *Y; int num_windows;
    int Hs; int out_L; fftwf_plan p_r2c; fftwf_plan p_c2r;

    float time_stretch = 2.0f;

    settup_vocoder(&time_buf, &win, &ifft_buf, &omega, &out, &norm, &new_data, &prev_phase, &sum_phase, 
                    &X, &Y, time_stretch, &num_windows, &Hs, &out_L, &p_r2c, &p_c2r);

    snd_pcm_sframes_t sent = 0;
    snd_pcm_sframes_t rcvd = 0;

    while (true)
    {
        // Capture PERIOD_FRAMES
        rcvd = 0;
        while (rcvd < PERIOD_FRAMES)
        {
            snd_pcm_sframes_t r = snd_pcm_readi(
                capture_handle,
                buffer + rcvd * CHANNELS,
                PERIOD_FRAMES - rcvd);
            if (r < 0)
            {
                r = xrun_recover(capture_handle, (int)r);
                if (r < 0)
                {
                    cerr << "capture recover failed\n";
                    goto out;
                }
                continue;
            }
            rcvd += r;
        }

        /* Run phase vo */
        phase_vocoder(buffer, time_buf, win, ifft_buf, omega, out, norm, new_data, prev_phase, sum_phase, X, Y, 
                      num_windows, Hs, out_L, p_r2c, p_c2r);

        // Playback PERIOD_FRAMES
        sent = 0;
        while (sent < PERIOD_FRAMES)
        {
            snd_pcm_sframes_t w = snd_pcm_writei(
                playback_handle,
                new_data + sent * CHANNELS,
                PERIOD_FRAMES - sent);
            if (w < 0)
            {
                w = xrun_recover(playback_handle, (int)w);
                if (w < 0)
                {
                    cerr << "playback recover failed\n";
                    goto out;
                }
                continue;
            }
            sent += w;
        }
    }

out:
    snd_pcm_drop(playback_handle);
    snd_pcm_close(playback_handle);
    snd_pcm_close(capture_handle);
    return 0;
}
