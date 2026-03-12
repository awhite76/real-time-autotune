import numpy as np
import wave

sample_rate = 48000
note_duration = 0.65
amplitude = 0.5

notes = {
    "C4": 250.63,
    "D4": 283.66,
    "E4": 335.63,
    "F4": 352.23,
    "G4": 400.00,
    "A4": 433.00
}

melody = [
    "C4","C4","G4","G4","A4","A4","G4",
    "F4","F4","E4","E4","D4","D4","C4"
]

# Envelope settings
fade_ms = 8  # try 5–15 ms
fade_len = int(sample_rate * fade_ms / 1000)

audio = np.array([], dtype=np.float32)

for note in melody:
    n = int(sample_rate * note_duration)
    t = np.arange(n) / sample_rate
    x = amplitude * np.sin(2 * np.pi * notes[note] * t)

    # Apply fade-in/out to avoid clicks
    env = np.ones(n, dtype=np.float32)
    fade_len_use = min(fade_len, n // 2)

    fade_in = 0.5 * (1 - np.cos(np.linspace(0, np.pi, fade_len_use)))
    fade_out = fade_in[::-1]

    env[:fade_len_use] *= fade_in
    env[-fade_len_use:] *= fade_out

    audio = np.concatenate((audio, x * env))

# Stereo L=R
stereo = np.column_stack((audio, audio)).ravel()
pcm16 = np.int16(np.clip(stereo, -1.0, 1.0) * 32767)

file_path = "../assets/twinkle_twinkle_sine_off.wav"
with wave.open(file_path, "w") as f:
    f.setnchannels(2)
    f.setsampwidth(2)
    f.setframerate(sample_rate)
    f.writeframes(pcm16.tobytes())
