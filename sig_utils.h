#pragma once

#include <fftw3.h>
#include <complex>

#define REAL 0
#define IMAG 1

#define A0 0.35875
#define A1 0.48829
#define A2 0.14128
#define A3 0.01168

enum fft_type {
    FFT = -1,
    IFFT = 1,
};

double blk_harris(int i, int n);
double cpx_abs(fftw_complex* x, int i);

void fill_rand(std::complex<double>* x_n, int N);
