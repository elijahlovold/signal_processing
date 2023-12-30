#include <iostream>
#include <cmath>
#include <fstream>
#include <vector> 
#include <math.h>
#include <complex.h>
#include <cstdlib>
#include <ctime>

// blackmann harris consts
#define A0 0.35875
#define A1 0.48829
#define A2 0.14128
#define A3 0.01168

#define REAL 0
#define IMAG 1

#define SAMPLE_RATE 800.0

enum fft_type {
    FFT = -1,
    IFFT = 1,
};

double blk_harris(int i, int N);
void custom_DFT(std::complex<double>* x_n, std::complex<double>* X_k, int N, fft_type s = FFT);

void fft(const std::complex<double>* x_n, std::complex<double>* X_k, int N);

