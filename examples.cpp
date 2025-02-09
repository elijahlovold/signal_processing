#include <matplotlibcpp.h>
#include <fftw3.h>
#include <math.h>
#include <vector>
#include <sig_utils.h>
#include <transfer.h>

#define SAMPLE_RATE 800.0

using namespace std;
namespace plt = matplotlibcpp;
namespace tf = transfer_functions;

int main() {
    // auto ab_pair = tf::butter(6, 5);
    // std::vector<double> a = ab_pair.first;
    // std::vector<double> b = ab_pair.second;

    // tf::plot_ht(a, b, 40);


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

