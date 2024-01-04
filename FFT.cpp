#include <matplotlibcpp.h>
#include <iostream>
#include <fftw3.h>
#include <math.h>
#include <vector>
#include <complex>

// blackmann harris consts
#define A0 0.35875
#define A1 0.48829
#define A2 0.14128
#define A3 0.01168
 
#define REAL 0
#define IMAG 1

#define SAMPLE_RATE 800.0

double cpx_abs(fftw_complex* x, int i);
double blk_harris(int i, int n);

using namespace std;
namespace plt = matplotlibcpp;

double blk_harris(int i, int n);
double cpx_abs(fftw_complex* x, int i);

int main() {
    // Set the size of the input array
    const int n = SAMPLE_RATE*5;

    vector<double> x_plot(n);
    vector<double> y_plot(n);

    fftw_complex x[n];
    fftw_complex x_win[n];
    fftw_complex y[n];

    for (int i = 0; i < n; ++i) {
        x[i][REAL] = sin(2.0*M_PI*i/SAMPLE_RATE) + sin(100*2.0*M_PI*i/SAMPLE_RATE); 
        // x[i][REAL] = sin(M_PI*(i-n/2)/SAMPLE_RATE)/(M_PI*(i-n/2)/SAMPLE_RATE); 
        x[i][IMAG] = 0;
    }
    
    for (int i = 0; i < n; ++i) {
        x_win[i][REAL] = x[i][REAL] * blk_harris(i, n); 
        x_win[i][IMAG] = 0;
    }

    fftw_plan plan = fftw_plan_dft_1d(n, x_win, y, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (int i = 0; i < n; ++i) {
        x_plot[i] = x[i][REAL];
        y_plot[i] = cpx_abs(y, i);
    }

    plt::figure();
    plt::plot(x_plot);
    
    plt::figure();
    plt::plot(y_plot);
    
    plt::show();

    return 0;


}

double blk_harris(int i, int n){
    return A0 
         - A1 * cos(2.0*M_PI*i / (n - 1))
         + A2 * cos(4.0*M_PI*i / (n - 1)) 
         - A3 * cos(6.0*M_PI*i / (n - 1));
}

double cpx_abs(fftw_complex* x, int i){
    return sqrt(x[i][REAL]*x[i][REAL] + x[i][IMAG]*x[i][IMAG]);
}
