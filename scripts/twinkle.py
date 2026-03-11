# Generate "Twinkle Twinkle Little Star" as pure sine waves (STEREO)

import numpy as np
import wave

# Audio settings
sample_rate = 48000  # Hz
note_duration = 0.65  # slightly slower tempo
amplitude = 0.5

# Frequencies (C major scale, 4th octave)
notes = {
    "C4": 261.63,
    "D4": 293.66,
    "E4": 329.63,
    "F4": 349.23,
    "G4": 392.00,
    "A4": 440.00
}

# Melody
melody = [
    "C4", "C4", "G4", "G4", "A4", "A4", "G4",
    "F4", "F4", "E4", "E4", "D4", "D4", "C4"
]

# Generate mono waveform first
audio = np.array([], dtype=np.float32)

for note in melody:
    t = np.linspace(0, note_duration, int(sample_rate * note_duration), endpoint=False)
    sine_wave = amplitude * np.sin(2 * np.pi * notes[note] * t)
    audio = np.concatenate((audio, sine_wave))

# Convert to stereo by duplicating channel (L = R)
stereo_audio = np.column_stack((audio, audio)).ravel()

# Convert to 16-bit PCM
audio_int16 = np.int16(stereo_audio * 32767)

# Save WAV file
file_path = "../assets/twinkle_twinkle_sine.wav"
with wave.open(file_path, 'w') as wav_file:
    wav_file.setnchannels(2)      # Stereo
    wav_file.setsampwidth(2)      # 16-bit
    wav_file.setframerate(sample_rate)
    wav_file.writeframes(audio_int16.tobytes())

file_path
