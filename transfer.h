#pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector> 
#include <math.h>
#include <complex.h>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

namespace transfer_functions {
    enum graph_type{
        ABS = true,
        LOG = false
    };
    enum freq_type{
        HZ = true,
        RAD = false
    };
    enum TF_type{
        ANA = true,
        DIG = false
    };
    enum filter_type{
        LPF,
        HPF, 
        BPF
    };

    void plot_ht(const std::vector<double>& b, const std::vector<double>& a, int N=50);
    void plot_TF(const std::vector<double>& b, const std::vector<double>& a, TF_type H_=DIG, freq_type f=HZ, graph_type t=ABS, int N=10000, double sample_rate=8000);
    void plot_pz(const std::vector<double>& b, const std::vector<double>& a, TF_type H_=DIG);
    void plot_t(const std::vector<double>& x_n, double sample_rate);
    void print_TF(const std::vector<double>& b, const std::vector<double>& a, TF_type H_=DIG);

    std::pair<std::vector<double>, std::vector<double>> butter(int N, double Omega_c=1, filter_type f=LPF);
    double butter_N(double WP, double Rp, double WS, double As, bool DB=false);

    std::pair<std::vector<double>, std::vector<double>> bilinear(const std::vector<double>& b_s, const std::vector<double>& a_s, double T=1);
    std::pair<std::vector<double>, std::vector<double>> zpk_to_ab(const std::vector<std::complex<double>>& z, const std::vector<std::complex<double>>& p, double k=1);
    std::vector<double> filt(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x);

    void print_poly(const std::vector<double>& b, const char var='x', bool inv=false);
    double poly_eval(const std::vector<double>& a, double x);
    std::complex<double> poly_eval(const std::vector<double>& a, std::complex<double> x);
    std::vector<double> poly_add(const std::vector<double>& b, const std::vector<double>& a);
    std::vector<double> poly_scale(const std::vector<double>& b, double x);
    std::vector<double> poly_multi(const std::vector<double>& b, const std::vector<double>& a);
    std::vector<double> poly_pow(const std::vector<double>& a, int n);

    std::vector<std::complex<double>> abrth_roots(const std::vector<double>& p, double tolerance=1e-8, int maxIterations=1e5);
}
