#include <iostream>
#include <fftw3.h>
#include <cmath>
#include <fstream>
#include <vector> 

// blackmann harris consts
#define A0 0.35875
#define A1 0.48829
#define A2 0.14128
#define A3 0.01168

#define REAL 0
#define IMAG 1

#define SAMPLE_RATE 8000.0

double cpx_abs(fftw_complex* x, int i);
double blk_harris(int i, int n);

void fft_arr_to_arr(fftw_complex* x, fftw_complex* y, int n);
void fft_arr_to_file(fftw_complex* x, int n);
bool fft_file_to_file(std::string filename, int n);
bool print_to_file(std::string filename, fftw_complex* x, int n);