#include <sig_utils.h>
#include <math.h>
#include <ctime>

double blk_harris(int i, int n){
    return A0 
         - A1 * cos(2.0*M_PI*i / (n - 1))
         + A2 * cos(4.0*M_PI*i / (n - 1)) 
         - A3 * cos(6.0*M_PI*i / (n - 1));
}

double cpx_abs(fftw_complex* x, int i){
    return sqrt(x[i][REAL]*x[i][REAL] + x[i][IMAG]*x[i][IMAG]);
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
