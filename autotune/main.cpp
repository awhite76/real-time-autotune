#include "main.hpp"
#include "pitch.hpp"
#include "time_stretch.hpp"
#include "pv.hpp"
#include "util.hpp"

int xrun_recover(snd_pcm_t *handle, int err)
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

bool set_hw_params(snd_pcm_t *handle, snd_pcm_stream_t stream)
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

// Compute frames for a T_MS window
static inline uint32_t framesForMs(uint32_t sampleRate, float ms)
{
    double frames = (double)sampleRate * ((double)ms / 1000.0);
    uint32_t out = (uint32_t)(frames + 0.5);
    return max<uint32_t>(out, 1);
}

PitchSeries buildPitchSeries_Tms(const StereoWavI16 &wav, float T_MS)
{
    if (wav.sampleRate == 0)
        throw runtime_error("Invalid sample rate");
    if (wav.left.size() != wav.right.size())
        throw runtime_error("L/R size mismatch");
    if (wav.left.empty())
        throw runtime_error("Empty WAV buffers");
    if (T_MS <= 0.0f)
        throw runtime_error("T_MS must be > 0");

    PitchSeries ps;
    ps.sampleRate = wav.sampleRate;
    ps.windowMs = T_MS;
    ps.windowFrames = framesForMs(wav.sampleRate, T_MS);

    const uint32_t N = ps.windowFrames;
    const size_t totalFrames = wav.left.size();
    const size_t chunks = totalFrames / N;

    ps.leftHz.resize(chunks);
    ps.rightHz.resize(chunks);
    ps.leftConf.resize(chunks);
    ps.rightConf.resize(chunks);

    // Yin detectors sized to the window length
    Yin yinL((int)N);
    Yin yinR((int)N);

    // Scratch buffers for each chunk (YIN expects contiguous int16_t*)
    vector<int16_t> bufL(N);
    vector<int16_t> bufR(N);

    for (size_t k = 0; k < chunks; ++k)
    {
        const size_t offset = k * (size_t)N;

        // copy chunk into contiguous buffers
        copy_n(wav.left.begin() + offset, N, bufL.begin());
        copy_n(wav.right.begin() + offset, N, bufR.begin());

        const float f0L = yinL.getPitch(bufL.data());
        const float cL = yinL.getProbability();

        const float f0R = yinR.getPitch(bufR.data());
        const float cR = yinR.getProbability();

        ps.leftHz[k] = f0L;
        ps.leftConf[k] = cL;
        ps.rightHz[k] = f0R;
        ps.rightConf[k] = cR;
    }

    return ps;
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

    /****************** Load .wav *********************/
    StereoWavI16 wav = loadStereoWav_i16("../assets/440Hz.wav");

    // If you want to ensure it matches your pipeline rate:
    if (wav.sampleRate != SAMPLE_RATE)
    {
        cerr << "Warning: input.wav sample rate is " << wav.sampleRate
             << " but pipeline expects " << SAMPLE_RATE << "\n";
    }

    // Build pitch time series (per T_MS chunk)
    PitchSeries pitchTS = buildPitchSeries_Tms(wav, (float)T_MS);

    cerr << "Loaded input.wav: frames=" << wav.left.size()
         << " windowFrames=" << pitchTS.windowFrames
         << " chunks=" << pitchTS.size() << "\n";

    // Example quick access:
    if (pitchTS.size() > 0)
    {
        size_t idx = pitchTS.indexForMs(250.0f); // chunk at ~250ms
        cerr << "Pitch @250ms: L=" << pitchTS.leftHz[idx] << "Hz"
             << " R=" << pitchTS.rightHz[idx] << "Hz\n";
    }

    /****************** audio buffers *****************/
    int16_t buffer[BUFFER_FRAMES];
    memset(buffer, 0, sizeof(buffer));

    // Per-channel mono buffers (PERIOD_FRAMES samples each)
    int16_t left[PERIOD_FRAMES];
    int16_t right[PERIOD_FRAMES];
    memset(left, 0, sizeof(left));
    memset(right, 0, sizeof(right));

    // Two independent YIN detectors (each holds its own yinBuffer/probability)
    Yin yinL(PERIOD_FRAMES);
    Yin yinR(PERIOD_FRAMES);

    /************* Time stretch config *************/
    TimeStretchResampler rs;
    time_stretch_init(rs, SAMPLE_RATE, 5);
    static int16_t rs_out[BUFFER_FRAMES];

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

    cout << "Looping " << T_MS << "ms chunks: "
         << PERIOD_FRAMES << " frames @ " << SAMPLE_RATE << " Hz, "
         << CHANNELS << " ch\n"
         << "Capture: " << cap_dev << "  Playback: " << pb_dev << "\n";

    /* Phase Vo settup */

    float *time_buf;
    float *win;
    float *ifft_buf;
    float *omega;
    float *out;
    float *norm;
    int16_t *new_data;
    float *prev_phase;
    float *sum_phase;
    fftwf_complex *X;
    fftwf_complex *Y;
    int num_windows;
    int Hs;
    int out_L;
    fftwf_plan p_r2c;
    fftwf_plan p_c2r;

    // 0.84 to 1.19 ~ +-3 semitones
    float time_stretch = 0.75f;

    cout << "Prevocoder\n";

    int vor = settup_vocoder(&time_buf, &win, &ifft_buf, &omega, &out, &norm, &new_data, &prev_phase, &sum_phase,
                             &X, &Y, time_stretch, &num_windows, &Hs, &out_L, &p_r2c, &p_c2r);

    printf("Vocoder settup: %d\n", vor);

    cout << "Post Vocoder\n";

    snd_pcm_sframes_t sent = 0;
    snd_pcm_sframes_t rcvd = 0;

    while (true)
    {
        // Capture PERIOD_FRAMES
        rcvd = 0;
        while (rcvd < PERIOD_FRAMES)
        {
            cout << "In the reading portion\n";

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

        /**************** Yin pitch detection ***************/
        // deinterleave_stereo_i16(buffer, left, right, PERIOD_FRAMES);

        // float f0L = yinL.getPitch(left);
        // float cL = yinL.getProbability();

        // float f0R = yinR.getPitch(right);
        // float cR = yinR.getProbability();

        // float f0Best = (cL >= cR) ? f0L : f0R;
        // float cBest = (cL >= cR) ? cL : cR;
        // const char *chBest = (cL >= cR) ? "L" : "R";

        // static int printCountdown = 0;
        // if (++printCountdown >= 10)
        // {
        //     printCountdown = 0;
        //     if (f0Best > 0.0f)
        //         cerr << "best(" << chBest << "): f0=" << f0Best << " Hz conf=" << cBest << "\n";
        //     else
        //         cerr << "best(" << chBest << "): f0=none conf=" << cBest << "\n";
        // }

        /* Run phase vo */
        memset(out, 0, (size_t)out_L * NUM_CHANNELS * sizeof(float));
        memset(norm, 0, (size_t)out_L * sizeof(float));
        phase_vocoder(buffer, time_buf, win, ifft_buf, omega, out, norm, new_data, prev_phase, sum_phase, X, Y,
                      num_windows, Hs, out_L, p_r2c, p_c2r);

        int outFrames = time_stretch_process(
            rs,
            new_data,
            out_L, // frames per channel available in new_data
            rs_out,
            PERIOD_FRAMES, // we want exactly one period for ALSA
            time_stretch); // ratio

        if (outFrames == 0)
        {
            cerr << "Speex resample failed\n";
            // fallback: play something sane
            memcpy(rs_out, buffer, sizeof(rs_out));
            outFrames = PERIOD_FRAMES;
        }

        if (outFrames < PERIOD_FRAMES)
        {
            memset(rs_out + outFrames * CHANNELS, 0,
                   (PERIOD_FRAMES - outFrames) * CHANNELS * sizeof(int16_t));
            outFrames = PERIOD_FRAMES;
        }

        deinterleave_stereo_i16(rs_out, left, right, PERIOD_FRAMES);

        float f0L = yinL.getPitch(left);
        float cL = yinL.getProbability();

        float f0R = yinR.getPitch(right);
        float cR = yinR.getProbability();

        float f0Best = (cL >= cR) ? f0L : f0R;
        float cBest = (cL >= cR) ? cL : cR;
        const char *chBest = (cL >= cR) ? "L" : "R";

        static int printCountdown = 0;
        if (++printCountdown >= 10)
        {
            printCountdown = 0;
            if (f0Best > 0.0f)
                cerr << "best(" << chBest << "): f0=" << f0Best << " Hz conf=" << cBest << "\n";
            else
                cerr << "best(" << chBest << "): f0=none conf=" << cBest << "\n";
        }

        // Playback PERIOD_FRAMES
        sent = 0;
        while (sent < PERIOD_FRAMES)
        {

            cout << "In the writing portion\n";
            snd_pcm_sframes_t w = snd_pcm_writei(
                playback_handle,
                rs_out + sent * CHANNELS,
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
    time_stretch_destroy(rs);
    return 0;
}
