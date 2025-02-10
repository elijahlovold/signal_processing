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
    // Set the size of the input array
    const int n = SAMPLE_RATE*5;

    vector<double> x_plot(n);
    vector<double> y_plot(n);

    // fftw_complex x[n];
    // fftw_complex x_win[n];
    // fftw_complex y[n];

    vector<double> x(n);
    vector<double> y(n);

    for (int i = 0; i < n; ++i) {
        x[i] = sin(2.0*M_PI*i/SAMPLE_RATE) + sin(100*2.0*M_PI*i/SAMPLE_RATE); 
        // x[i][REAL] = sin(2.0*M_PI*i/SAMPLE_RATE) + sin(100*2.0*M_PI*i/SAMPLE_RATE); 
        // x[i][IMAG] = 0;
    }

    std::vector<double> a = {1, 2, 1};
    std::vector<double> b = {684.8, -1293.6, 612.8};

    tf::print_TF(a, b);
    tf::plot_TF(a, b);
    tf::plot_ht(a, b, 500);

    y = tf::filt(a, b, x);

    int N = x.size();

    std::vector<double> ind(N);

    for (double i = 0; i < N; i++){
        ind[i] = i/SAMPLE_RATE;
    }

    plt::figure();
    plt::plot(ind, x);
    plt::plot(ind, y);
    plt::title("t-sig");
    plt::xlabel("Time (s)");
    plt::ylabel("Magnitude");
    plt::grid(true);

    plt::show();

    // for (int i = 0; i < n; ++i) {
    //     x_win[i][REAL] = x[i][REAL] * blk_harris(i, n); 
    //     x_win[i][IMAG] = 0;
    // }

    // fftw_plan plan = fftw_plan_dft_1d(n, x_win, y, FFTW_FORWARD, FFTW_ESTIMATE);
    // fftw_execute(plan);
    // fftw_destroy_plan(plan);

    // for (int i = 0; i < n; ++i) {
    //     x_plot[i] = x[i][REAL];
    //     y_plot[i] = cpx_abs(y, i);
    // }

    // plt::figure();
    // plt::plot(x_plot);

    // plt::figure();
    // plt::plot(y_plot);

    // plt::show();

    // return 0;


}

