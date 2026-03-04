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

#include "util.hpp"         // loadStereoWav_i16, interleaveStereo, writeStereoWav_i16_interleaved
#include "time_stretch.hpp" // TimeStretchResampler + time_stretch_*
#include "pv.hpp"           // settup_vocoder, phase_vocoder, NUM_CHANNELS, NUM_FRAMES

using namespace std;

static string make_outname(const string &inPath, const string &tag, float s)
{
    auto dot = inPath.find_last_of('.');
    string stem = (dot == string::npos) ? inPath : inPath.substr(0, dot);
    ostringstream ss;
    ss << stem << "_" << tag << "_s" << fixed << setprecision(2) << s;
    string name = ss.str();
    for (char &c : name)
        if (c == '.')
            c = 'p';
    name += ".wav";
    return name;
}

int main(int argc, char **argv)
{
    try
    {
        string inPath = (argc > 1) ? argv[1] : "../assets/440Hz.wav";

        // Ensure output directory exists
        const char *outdir = "../out/pv_out";
        mkdir(outdir, 0755); // ignore errors if exists

        // Load file
        StereoWavI16 wav = loadStereoWav_i16(inPath);
        if (wav.sampleRate == 0)
            throw runtime_error("Invalid WAV sample rate");
        if (wav.left.empty())
            throw runtime_error("Empty WAV");

        // Interleave LRLR...
        vector<int16_t> inter = interleaveStereo(wav);
        const int totalFrames = static_cast<int>(wav.left.size());

        // We'll process input in chunks of NUM_FRAMES frames (vocoder expects that)
        const int chunkFrames = NUM_FRAMES; // from pv.hpp
        const size_t samplesPerChunk = (size_t)chunkFrames * (size_t)NUM_CHANNELS;

        // Setup vocoder for chunk-sized processing (allocate buffers sized for chunkFrames)
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

        // If your settup_vocoder takes input_frames, pass chunkFrames to size buffers appropriately.
        // If your version does not take input_frames, call it as you normally do.
        int rc = settup_vocoder(&time_buf, &win, &ifft_buf, &omega, &out, &norm, &new_data, &prev_phase, &sum_phase,
                                &X, &Y, chunkFrames, &num_windows, &Hs_alloc, &max_out_L, &p_r2c, &p_c2r);
        if (rc != 0)
            throw runtime_error("settup_vocoder failed");

        cerr << "Vocoder setup: chunkFrames=" << chunkFrames << " num_windows=" << num_windows
             << " alloc_out_L=" << max_out_L << "\n";

        // Init speex resampler
        TimeStretchResampler rs;
        if (!time_stretch_init(rs, wav.sampleRate, /*quality=*/5))
        {
            cerr << "Warning: time_stretch_init failed (resampling disabled)\n";
        }

        // Buffer we pass to phase_vocoder for each chunk (interleaved samples)
        vector<int16_t> chunkBuf(samplesPerChunk);

        // For each stretch factor s, run vocoder repeatedly over input chunks and concatenate outputs
        for (float s = 0.40f; s <= 2.5001f; s += 0.30f)
        {
            cout << "=== Processing s=" << fixed << setprecision(2) << s << " ===\n";

            // clear vocoder state so each s run starts fresh (optional)
            if (prev_phase)
                memset(prev_phase, 0, (size_t)NUM_CHANNELS * (size_t)FREQ_BINS * sizeof(float));
            if (sum_phase)
                memset(sum_phase, 0, (size_t)NUM_CHANNELS * (size_t)FREQ_BINS * sizeof(float));

            // storage for concatenated phase-vocoder output (interleaved int16)
            vector<int16_t> pv_concat;
            pv_concat.reserve((size_t)totalFrames * (size_t)NUM_CHANNELS / 2); // heuristic reserve

            int readPos = 0;
            while (readPos < totalFrames)
            {
                // copy next chunkFrames frames (zero-pad tail)
                int avail = totalFrames - readPos;
                int toCopy = std::min(avail, chunkFrames);
                // interleaved index = readPos * NUM_CHANNELS
                memcpy(chunkBuf.data(), inter.data() + (size_t)readPos * NUM_CHANNELS, (size_t)toCopy * NUM_CHANNELS * sizeof(int16_t));
                if (toCopy < chunkFrames)
                {
                    // zero pad rest
                    memset(chunkBuf.data() + (size_t)toCopy * NUM_CHANNELS, 0, (size_t)(chunkFrames - toCopy) * NUM_CHANNELS * sizeof(int16_t));
                }

                // Zero out float accum buffers (out, norm) before call
                if (out && norm)
                {
                    memset(out, 0, (size_t)max_out_L * (size_t)NUM_CHANNELS * sizeof(float));
                    memset(norm, 0, (size_t)max_out_L * sizeof(float));
                }

                // Call phase vocoder on this chunk; it will update out_L for this chunk
                int out_L = 0;
                int pv_err = phase_vocoder(chunkBuf.data(), time_buf, win, ifft_buf, omega,
                                           out, norm, new_data, prev_phase, sum_phase, X, Y,
                                           s, &out_L, num_windows, p_r2c, p_c2r);
                if (pv_err != 0)
                {
                    cerr << "phase_vocoder failed for s=" << s << " (err=" << pv_err << ")\n";
                    break;
                }

                // Append new_data[0 .. out_L*NUM_CHANNELS-1] to pv_concat
                size_t samplesProduced = (size_t)out_L * (size_t)NUM_CHANNELS;
                pv_concat.insert(pv_concat.end(), new_data, new_data + samplesProduced);

                // Advance read position by chunkFrames (non-overlapping chunks).
                // If you prefer sliding processing, advance by ANALYSIS_HOP or another hop.
                readPos += chunkFrames;
            } // end per-chunk loop

            // Write concatenated pv result
            string pv_outname = string(outdir) + "/pv_out_s" + to_string((int)lroundf(s * 100.0f)) + ".wav";
            writeStereoWav_i16_interleaved(pv_outname, wav.sampleRate, pv_concat.data(), (uint32_t)(pv_concat.size() / NUM_CHANNELS));
            cout << "Wrote " << pv_outname << " frames=" << (pv_concat.size() / NUM_CHANNELS) << "\n";

            // Now resample the entire concatenated pv output by 's' to produce pitched result
            if (!rs.st)
            {
                cerr << "Resampler not initialized; skipping resample for s=" << s << "\n";
                continue;
            }

            // conservative output capacity estimate
            int inFrames_pv = (int)(pv_concat.size() / NUM_CHANNELS);
            int outCapacityFrames = (int)ceil((double)inFrames_pv * max(1.0, (double)s)) + 1024;
            vector<int16_t> rs_out((size_t)outCapacityFrames * NUM_CHANNELS);

            int outFrames = time_stretch_process(rs, pv_concat.data(), inFrames_pv, rs_out.data(), outCapacityFrames, s);
            if (outFrames <= 0)
            {
                cerr << "time_stretch_process failed for s=" << s << "\n";
                continue;
            }

            string pitch_outname = string(outdir) + "/pv_pitch_s" + to_string((int)lroundf(s * 100.0f)) + ".wav";
            writeStereoWav_i16_interleaved(pitch_outname, wav.sampleRate, rs_out.data(), (uint32_t)outFrames);
            cout << "Wrote pitched file " << pitch_outname << " frames=" << outFrames << "\n";
        } // end s sweep

        // cleanup
        time_stretch_destroy(rs);

        // free allocated stuff (mirror settup_vocoder allocations)
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