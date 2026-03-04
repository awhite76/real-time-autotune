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

#include "util.hpp"         // loadStereoWav_i16
#include "time_stretch.hpp" // TimeStretchResampler + time_stretch_*

void test_resample(const std::string &inputWavPath)
{
    StereoWavI16 wav = loadStereoWav_i16(inputWavPath);
    if (wav.sampleRate == 0)
        throw std::runtime_error("Invalid WAV sample rate");
    if (wav.left.empty())
        throw std::runtime_error("Empty WAV");

    // Interleaved input (LRLR...)
    std::vector<int16_t> inInter = interleaveStereo(wav);
    const int inFrames = static_cast<int>(wav.left.size());

    TimeStretchResampler r{};
    if (!time_stretch_init(r, wav.sampleRate, /*quality=*/5))
    {
        throw std::runtime_error("time_stretch_init failed");
    }

    // Sweep s from 0.40 to 2.5 in steps of 0.3
    for (float s = 0.40f; s <= 2.5001f; s += 0.30f)
    {
        // Worst-case output frames: if s is small, output grows roughly like 1/s.
        // Allocate a conservative bound (+a little headroom).
        const int outCapacityFrames = static_cast<int>(std::ceil(inFrames / std::max(s, 0.01f))) + 1024;

        std::vector<int16_t> outInter(static_cast<size_t>(outCapacityFrames) * 2);

        int outFrames = time_stretch_process(
            r,
            inInter.data(),
            inFrames, // frames per channel
            outInter.data(),
            outCapacityFrames, // frames per channel capacity
            s);

        if (outFrames <= 0)
        {
            std::cerr << "Resample failed at s=" << s << "\n";
            continue;
        }

        std::string outName = "../out/resampled_" + std::to_string(static_cast<int>(s * 100)) + "pct.wav";
        writeStereoWav_i16_interleaved(
            outName,
            wav.sampleRate,
            outInter.data(),
            static_cast<uint32_t>(outFrames));

        std::cout << "Wrote " << outName
                  << "  inFrames=" << inFrames
                  << "  outFrames=" << outFrames
                  << "  s=" << std::fixed << std::setprecision(2) << s
                  << "\n";
    }

    time_stretch_destroy(r);
}

int main(int argc, char **argv)
{
    try
    {
        std::string in = (argc > 1) ? argv[1] : "../assets/440Hz.wav";
        test_resample(in);
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_resample failed: " << e.what() << "\n";
        return 1;
    }
}