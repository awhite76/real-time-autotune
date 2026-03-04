// test_pv.cpp
#include <cstdint>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h> // mkdir
#include <sys/types.h>

#include "util.hpp"
#include "time_stretch.hpp"
#include "pv.hpp"

using namespace std;

int main(int argc, char **argv)
{
    try
    {
        string inPath = (argc > 1) ? argv[1] : "../assets/440Hz.wav";

        const char *outdir = "../out/pv_out";
        mkdir(outdir, 0755);

        StereoWavI16 wav = loadStereoWav_i16(inPath);
        if (wav.sampleRate == 0)
            throw runtime_error("Invalid WAV sample rate");
        if (wav.left.empty())
            throw runtime_error("Empty WAV");

        vector<int16_t> inter = interleaveStereo(wav);
        const int totalFrames = (int)wav.left.size();

        // Chunk size the PV code expects (pv.cpp uses NUM_FRAMES internally)
        const int chunkFrames = NUM_FRAMES;
        const size_t samplesPerChunk = (size_t)chunkFrames * (size_t)NUM_CHANNELS;

        // --- Vocoder setup (signature per your pv.hpp) ---
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
        int Hs_alloc = 0;
        int max_out_L = 0;
        fftwf_plan p_r2c = nullptr;
        fftwf_plan p_c2r = nullptr;

        // IMPORTANT: no &chunkFrames here; your current API does not accept it.
        int rc = settup_vocoder(&time_buf, &win, &ifft_buf, &omega,
                                &out, &norm, &new_data,
                                &prev_phase, &sum_phase,
                                &X, &Y,
                                &num_windows, &Hs_alloc,
                                &max_out_L, &p_r2c, &p_c2r);

        if (rc != 0)
            throw runtime_error("settup_vocoder failed");

        cerr << "Vocoder setup ok: chunkFrames(NUM_FRAMES)=" << chunkFrames
             << " num_windows=" << num_windows
             << " max_out_L=" << max_out_L << "\n";

        // --- Resampler ---
        TimeStretchResampler rs{};
        if (!time_stretch_init(rs, wav.sampleRate, 5))
            cerr << "Warning: time_stretch_init failed (resampling disabled)\n";

        vector<int16_t> chunkBuf(samplesPerChunk);

        for (float s = 0.40f; s <= 2.5001f; s += 0.30f)
        {
            cout << "=== s=" << fixed << setprecision(2) << s << " ===\n";

            // Reset PV phase state per run so each output file is independent
            if (prev_phase)
                memset(prev_phase, 0, (size_t)NUM_CHANNELS * (size_t)FREQ_BINS * sizeof(float));
            if (sum_phase)
                memset(sum_phase, 0, (size_t)NUM_CHANNELS * (size_t)FREQ_BINS * sizeof(float));

            vector<int16_t> pv_concat;
            pv_concat.reserve((size_t)totalFrames * (size_t)NUM_CHANNELS); // conservative

            int readPos = 0;
            int chunkCount = 0;

            while (readPos < totalFrames)
            {
                int avail = totalFrames - readPos;
                int toCopy = std::min(avail, chunkFrames);

                // Copy interleaved frames into chunkBuf
                memcpy(chunkBuf.data(),
                       inter.data() + (size_t)readPos * NUM_CHANNELS,
                       (size_t)toCopy * NUM_CHANNELS * sizeof(int16_t));

                // Zero-pad tail
                if (toCopy < chunkFrames)
                {
                    memset(chunkBuf.data() + (size_t)toCopy * NUM_CHANNELS,
                           0,
                           (size_t)(chunkFrames - toCopy) * NUM_CHANNELS * sizeof(int16_t));
                }

                // Clear accumulators for this chunk (BUT keep prev_phase/sum_phase!)
                memset(out, 0, (size_t)max_out_L * (size_t)NUM_CHANNELS * sizeof(float));
                memset(norm, 0, (size_t)max_out_L * sizeof(float));

                int out_L = 0;
                int pv_err = phase_vocoder(chunkBuf.data(), time_buf, win, ifft_buf, omega,
                                           out, norm, new_data,
                                           prev_phase, sum_phase, X, Y,
                                           s, &out_L, num_windows, p_r2c, p_c2r);
                if (pv_err != 0)
                    throw runtime_error("phase_vocoder failed");

                // Append produced samples
                pv_concat.insert(pv_concat.end(),
                                 new_data,
                                 new_data + (size_t)out_L * (size_t)NUM_CHANNELS);

                readPos += chunkFrames;
                chunkCount++;
            }

            const uint32_t pvFrames = (uint32_t)(pv_concat.size() / NUM_CHANNELS);
            string pv_outname = string(outdir) + "/pv_" + to_string((int)lroundf(s * 100.0f)) + "pct.wav";
            writeStereoWav_i16_interleaved(pv_outname, wav.sampleRate, pv_concat.data(), pvFrames);
            cout << "Wrote " << pv_outname << " frames=" << pvFrames
                 << " chunks=" << chunkCount << "\n";

            // Resample by same s to get pitch shift (per your pipeline)
            if (!rs.st)
            {
                cerr << "Resampler not initialized; skipping resample\n";
                continue;
            }

            const int inFrames_pv = (int)pvFrames;

            // output frames often scale ~inFrames*s (based on your earlier use)
            int outCapacityFrames = (int)ceil((double)inFrames_pv * std::max(1.0, (double)s)) + 1024;
            vector<int16_t> rs_out((size_t)outCapacityFrames * (size_t)NUM_CHANNELS);

            int outFrames = time_stretch_process(rs,
                                                 pv_concat.data(),
                                                 inFrames_pv,
                                                 rs_out.data(),
                                                 outCapacityFrames,
                                                 s);
            if (outFrames <= 0)
                throw runtime_error("time_stretch_process failed");

            string pitch_outname = string(outdir) + "/pv_resampled_" + to_string((int)lroundf(s * 100.0f)) + "pct.wav";
            writeStereoWav_i16_interleaved(pitch_outname, wav.sampleRate, rs_out.data(), (uint32_t)outFrames);
            cout << "Wrote " << pitch_outname << " frames=" << outFrames << "\n";
        }

        time_stretch_destroy(rs);

        // cleanup (mirror settup_vocoder allocations)
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

        cout << "Done.\n";
        return 0;
    }
    catch (const exception &e)
    {
        cerr << "test_pv failed: " << e.what() << "\n";
        return 1;
    }
}