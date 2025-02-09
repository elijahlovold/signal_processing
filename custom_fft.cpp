#include "custom_fft.h"
#include <sig_utils.h>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;

#ifndef IS_MAIN_FILE
int main() {
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Set the size of the input array
    const int N = SAMPLE_RATE*5;
    // const int N = 2 << 10;

    std::complex<double> X_k[N];
    
    std::complex<double> x_n[N];
    // for (int i = 0; i < N; ++i) {
    //     x_n[i].real(sin(2.0*M_PI*i/SAMPLE_RATE) + sin(100*2.0*M_PI*i/SAMPLE_RATE));
    //     x_n[i].imag(0);
    // }
    
    fill_rand(x_n, N);

    std::complex<double> x_n_win[N];
    for (int i = 0; i < N; ++i) {
        x_n_win[i].real(x_n[i].real() * blk_harris(i, N));
        x_n_win[i].imag(0);
    }

    // custom_DFT(x_n_win, X_k, N);        // DFT
    custom_DFT(x_n, X_k, N);        // DFT

    // add harmonics 
    X_k[int(60.0*N/SAMPLE_RATE)].real(10000);
    X_k[int(120.0*N/SAMPLE_RATE)].real(3000);
    X_k[int(180.0*N/SAMPLE_RATE)].real(1000);
    
    custom_DFT(X_k, x_n, N, IFFT);    // IDFT

    // fft(x_n, X_k, N);        // FFT

    std::vector<double> time_domain(N);
    std::vector<double> freq_domain(N);
    // std::vector<double> times = gen_times(SAMPLE_RATE, N);
    // std::vector<double> freqs = gen_freqs(SAMPLE_RATE, N);

    for (int i = 0; i < N; ++i) {
        time_domain[i] = x_n[i].real();
        freq_domain[i] = abs(X_k[i]);
    }


    plt::figure(1);
    plt::plot(gen_times(SAMPLE_RATE, N), time_domain);
    plt::title("t-sig");
    plt::xlabel("time (s)");
    plt::ylabel("Magnitude");
    plt::grid(true);

    plt::figure(2);
    plt::plot(gen_freqs(SAMPLE_RATE, N), freq_domain);
    plt::title("fft");
    plt::xlabel("Frequency (Hz)");
    plt::ylabel("Magnitude");
    
    plt::show();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("CPU Time Used: %f seconds\n", cpu_time_used);

    return 0;
}
#endif

void custom_DFT(std::complex<double>* x1, std::complex<double>* x2, int N, fft_type s){
    for (int k = 0; k < N; k++){
        x2[k] = std::complex<double> (0,0);
        for (int n = 0; n < N; n++){
            x2[k] += x1[n]*(std::polar(1.0, s*(2.0*M_PI*k*n)/N));
        }
        if (s == IFFT){   // scale if IDFT
            x2[k] *= 1.0/N;
        }
    }
}

void fft(const std::complex<double>* x_n, std::complex<double>* X_k, int N) {
    if (N <= 1) {
        X_k[0] = x_n[0];
    } else {
        // Divide step
        std::complex<double>* x_even = new std::complex<double>[N / 2];
        std::complex<double>* x_odd = new std::complex<double>[N / 2];

        for (int i = 0; i < N / 2; ++i) {
            x_even[i] = x_n[2 * i];
            x_odd[i] = x_n[2 * i + 1];
        }

        // Recursion step
        std::complex<double>* X_even = new std::complex<double>[N / 2];
        std::complex<double>* X_odd = new std::complex<double>[N / 2];

        fft(x_even, X_even, N / 2);
        fft(x_odd, X_odd, N / 2);

        // Combine step
        for (int k = 0; k < N / 2; ++k) {
            std::complex<double> t = std::polar(1.0, -2.0 * M_PI * k / N) * X_odd[k];
            X_k[k] = X_even[k] + t;
            X_k[k + N / 2] = X_even[k] - t;
        }

        // Clean up
        delete[] x_even;
        delete[] x_odd;
        delete[] X_even;
        delete[] X_odd;
    }
}



void fill_rand(std::complex<double>* x_n, int N){    
    std::srand(std::time(0));
    double position = double(std::rand())/RAND_MAX;  // Initial position
    for (int i = 0; i < N; ++i) {
        position += (double(std::rand())/RAND_MAX - 0.5);
        x_n[i].real(position);
        x_n[i].imag(0);
    }
}

std::vector<double> gen_freqs(double F_S, int N){
    std::vector<double> result(N);

    for(int i = 0; i < N; i++){
        result[i] = i*F_S/N;
    }

    return result;
}

std::vector<double> gen_times(double F_S, int N){
    std::vector<double> result(N);

    for(int i = 0; i < N; i++){
        result[i] = i/F_S;
    }

    return result;
}
