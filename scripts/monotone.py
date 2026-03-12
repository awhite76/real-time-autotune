import numpy as np
import wave

# Audio settings
sample_rate = 48000
note_duration = 0.65  # same tempo as before
amplitude = 0.5

# Single frequency
# 300 Hz 
single_freq = 300.0

melody_len = 14  # same beat pattern length

audio = np.array([], dtype=np.float32)

for _ in range(melody_len):
    t = np.linspace(0, note_duration, int(sample_rate * note_duration), endpoint=False)
    sine_wave = amplitude * np.sin(2 * np.pi * single_freq * t)
    audio = np.concatenate((audio, sine_wave))

# Stereo (L = R)
stereo_audio = np.column_stack((audio, audio)).ravel()
audio_int16 = np.int16(stereo_audio * 32767)

file_path = "../assets/twinkle_rhythm_one_note_300hz_stereo.wav"
with wave.open(file_path, "w") as wav_file:
    wav_file.setnchannels(2)
    wav_file.setsampwidth(2)
    wav_file.setframerate(sample_rate)
    wav_file.writeframes(audio_int16.tobytes())

file_path