#pragma once

#include <cmath>
#include <vector> 
#include <math.h>
#include <complex.h>
#include <cstdlib>
#include <ctime>
#include <sig_utils.h>

#define SAMPLE_RATE 800.0

void custom_DFT(std::complex<double>* x_n, std::complex<double>* X_k, int N, fft_type s = FFT);

void fft(const std::complex<double>* x_n, std::complex<double>* X_k, int N);

std::vector<double> gen_freqs(double F_S, int N);
std::vector<double> gen_times(double F_S, int N);
