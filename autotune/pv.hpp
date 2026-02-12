/* Functions to support phase vocoder implementation */

typedef struct {
    int N;              // FFT size / window length
    int H;              // hop size
    int num_frames;     // number of time frames
    float *window;      // length N
    float *time_buf;    // length N (input frame buffer)
    fftwf_complex *spec;// length num_frames * (N/2+1)
    fftwf_plan plan;    // r2c plan on time_buf -> temp_out
    fftwf_complex *temp_out; // length (N/2+1)
} STFT;

float win[]
float *time_buf, 
float *ifft_buf, 
fftwf_complex *X, 
fftwf_complex *Y, 
float *out, 
float *norm,
float *prev_phase, 
float *sum_phase
float *omega