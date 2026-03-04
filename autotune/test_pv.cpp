// test_pv.cpp
#include <cstdint>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h> // mkdir
#include <sys/types.h>

#include "util.hpp" // loadStereoWav_i16, interleaveStereo, writeStereoWav_i16_interleaved
#include "time_stretch.hpp"
#include "pv.hpp" // settup_vocoder, phase_vocoder, constants like NUM_CHANNELS, FREQ_BINS, WINDOW_SIZE

int main(int argc, char **argv)
{
    try
    {
        std::string inPath = (argc > 1) ? argv[1] : "../assets/440Hz.wav";

        // Ensure output directory exists
        const char *outdir = "../out/pv_out";
        mkdir(outdir, 0755); // ignore errors if exists

        // Load file
        StereoWavI16 wav = loadStereoWav_i16(inPath);
        if (wav.sampleRate == 0)
            throw std::runtime_error("Invalid WAV sample rate");
        if (wav.left.empty())
            throw std::runtime_error("Empty WAV");

        // Interleave LRLR...
        std::vector<int16_t> inter = interleaveStereo(wav);
        const int inFrames = static_cast<int>(wav.left.size());

        // Setup vocoder (allocs/win/fftw plans/etc)
        float *time_buf = nullptr;
        float *win = nullptr;
        float *ifft_buf = nullptr;
        float *omega = nullptr;
        float *out = nullptr;
        float *norm = nullptr;
        int16_t *new_data = nullptr;
        float *prev_phase = nullptr;
        float *sum_phase = nullptr;
        fftwf_complex *X = nullptr;
        fftwf_complex *Y = nullptr;
        int num_windows = 0;
        int Hs = 0;
        int out_L = 0;
        fftwf_plan p_r2c = nullptr;
        fftwf_plan p_c2r = nullptr;

        int rc = settup_vocoder(&time_buf, &win, &ifft_buf, &omega, &out, &norm, &new_data, &prev_phase, &sum_phase,
                                &X, &Y, &num_windows, &Hs, &out_L, &p_r2c, &p_c2r);
        if (rc != 0)
            throw std::runtime_error("settup_vocoder failed");

        std::cout << "Vocoder setup ok: num_windows=" << num_windows << "  initial_out_L=" << out_L << "\n";

        // Prepare Speex resampler instance (we will init/destroy per run via helper)
        TimeStretchResampler rs;
        if (!time_stretch_init(rs, wav.sampleRate, /*quality=*/5))
        {
            // Not fatal but warn
            std::cerr << "Warning: time_stretch_init failed (resampling disabled)\n";
        }

        // Sweep s
        for (float s = 0.40f; s <= 2.5001f; s += 0.30f)
        {
            std::cout << "=== Processing s=" << std::fixed << std::setprecision(2) << s << " ===\n";

            // Clear vocoder output buffers (out and norm and prev/sum phases)
            if (out && norm)
                memset(out, 0, (size_t)out_L * (size_t)NUM_CHANNELS * sizeof(float));
            if (norm)
                memset(norm, 0, (size_t)out_L * sizeof(float));
            if (prev_phase)
                memset(prev_phase, 0, (size_t)NUM_CHANNELS * (size_t)FREQ_BINS * sizeof(float));
            if (sum_phase)
                memset(sum_phase, 0, (size_t)NUM_CHANNELS * (size_t)FREQ_BINS * sizeof(float));

            // Call phase vocoder: produces new_data (int16 interleaved) and updates out_L based on s
            int pv_err = phase_vocoder(inter.data(), time_buf, win, ifft_buf, omega,
                                       out, norm, new_data, prev_phase, sum_phase, X, Y,
                                       s, &out_L, num_windows, p_r2c, p_c2r);
            if (pv_err != 0)
            {
                std::cerr << "phase_vocoder failed for s=" << s << " (err=" << pv_err << ")\n";
                continue;
            }

            // Write the raw phase-vocoder **time-stretched** output
            std::string pv_outname = std::string(outdir) + "/pv_" + std::to_string(static_cast<int>(s * 100)) + "pct.wav";
            writeStereoWav_i16_interleaved(pv_outname, wav.sampleRate, new_data, (uint32_t)out_L);
            std::cout << "Wrote phase-vocoder output: " << pv_outname << " frames=" << out_L << "\n";

            // Now resample the vocoder output by the same factor 's' to produce a pitch-shifted result.
            // If time_stretch_process maps inputFrames -> approx inputFrames * s outputFrames, allocate accordingly.
            if (!rs.st)
            {
                std::cerr << "Resampler not initialized; skipping resample for s=" << s << "\n";
                continue;
            }

            // estimate output capacity conservatively
            int outCapacityFrames = (int)std::ceil((double)out_L * std::max(1.0, (double)s)) + 1024;
            std::vector<int16_t> rs_out((size_t)outCapacityFrames * (size_t)NUM_CHANNELS);

            int outFrames = time_stretch_process(rs, new_data, out_L, rs_out.data(), outCapacityFrames, s);
            if (outFrames <= 0)
            {
                std::cerr << "time_stretch_process failed for s=" << s << "\n";
                continue;
            }

            std::string pitch_outname = std::string(outdir) + "/pv_resampled_" + std::to_string(static_cast<int>(s * 100)) + "pct.wav";
            writeStereoWav_i16_interleaved(pitch_outname, wav.sampleRate, rs_out.data(), (uint32_t)outFrames);
            std::cout << "Wrote pitched output: " << pitch_outname << " frames=" << outFrames << "\n";
        }

        // cleanup
        time_stretch_destroy(rs);

        // free/cleanup as settup_vocoder allocated (mirror its error path)
        if (new_data)
            free(new_data);
        if (omega)
            free(omega);
        if (prev_phase)
            free(prev_phase);
        if (sum_phase)
            free(sum_phase);
        if (out)
            free(out);
        if (norm)
            free(norm);
        if (p_r2c)
            fftwf_destroy_plan(p_r2c);
        if (p_c2r)
            fftwf_destroy_plan(p_c2r);
        if (win)
            fftwf_free(win);
        if (time_buf)
            fftwf_free(time_buf);
        if (ifft_buf)
            fftwf_free(ifft_buf);
        if (X)
            fftwf_free(X);
        if (Y)
            fftwf_free(Y);

        std::cout << "Done.\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_pv failed: " << e.what() << "\n";
        return 1;
    }
}