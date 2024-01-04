#include "transfer.h"

namespace transfer_functions {
    std::pair<std::vector<double>, std::vector<double>> bilinear(const std::vector<double>& b_s, const std::vector<double>& a_s, double T){
        int Nb = b_s.size() - 1;
        int Na = a_s.size() - 1;

        int NH = (Na > Nb ? Na : Nb);   // highest order 

        std::vector<double> b_z;    // num
        std::vector<double> a_z;    // den
        
        std::vector<double> bi_num = {2, -2};
        std::vector<double> bi_den = {T, T};

        std::vector<double> term;

        // compute new numenator
        for (int i = 0; i <= Nb; i++){
            term = poly_multi(poly_pow(bi_num, i), poly_pow(bi_den, (NH - i)));
            term = poly_scale(term, b_s[i]);
            b_z = poly_add(b_z, term);
        }

        // compute new denominator
        for (int i = 0; i <= Na; i++){
            term = poly_multi(poly_pow(bi_num, i), poly_pow(bi_den, (NH - i)));
            term = poly_scale(term, a_s[i]);
            a_z = poly_add(a_z, term);
        }

        // normalize
        double a_0 = a_z[0];
        for (int i = 0; i < b_z.size(); i++){
            b_z[i] /= a_0;
        }
        for (int i = 0; i < a_z.size(); i++){
            a_z[i] /= a_0;
        }

        print_TF(b_z, a_z);
        return std::make_pair(b_z, a_z);
    }

    std::complex<double> poly_eval(const std::vector<double>& a, std::complex<double> x){
        std::complex<double> y = 0;
        for (int n = 0; n < a.size(); n++){
            y += a[n]*pow(x, n);
        }
        return y;
    }

    double poly_eval(const std::vector<double>& a, double x){
        double y = 0;
        for (int n = 0; n < a.size(); n++){
            y += a[n]*pow(x, n);
        }
        return y;
    }

    std::vector<double> poly_add(const std::vector<double>& b, const std::vector<double>& a){
        int Na = a.size();
        int Nb = b.size();
        std::vector<double> z((Na > Nb ? Na : Nb));

        for (int i = 0; i < Na; i++){
            z[i] += a[i];
        }
        for (int i = 0; i < Nb; i++){
            z[i] += b[i];
        }
        return z;
    }

    std::vector<double> poly_scale(const std::vector<double>& b, double x){
        int N = b.size();
        std::vector<double> z(N);
        for (int i = 0; i < N; i++){
            z[i] = b[i]*x;
        }
        return z;
    }

    
    std::vector<double> poly_multi(const std::vector<double>& b, const std::vector<double>& a){
        int Nb = b.size();
        int Na = a.size();

        std::vector<double> z(Na + Nb - 1);

        for (int i = 0; i < Na; i++){
            for (int j = 0; j < Nb; j++){
                z[i + j] += a[i]*b[j];
            }
        }
        return z;
    }
     
    std::vector<double> poly_pow(const std::vector<double>& a, int n){
        std::vector<double> z;

        z = ((n > 0) ? a : std::vector<double>{1});
        for (int i = 0; i < n-1; i++){
            z = poly_multi(z, a);
        }
        return z;
    }

    std::vector<double> filt(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x){
        int N = x.size() - 1;
        std::vector<double> y(N, 0.0); 

        for (int i = 0; i < N; ++i){
            for (int n = 0; n < b.size(); ++n){
                y[i] += b[n] * (i >=n ? x[i-n] : 0);    // pad with zeros
            }

            for (int d = 1; d < a.size(); ++d){
                y[i] -= a[d] * (i >=d ? y[i-d] : 0);    // pad with zeros
            }
            y[i] /= a[0];
        }

        return y;
    }

    void plot_pz(const std::vector<double>& b, const std::vector<double>& a, TF_type H_){
        std::vector<std::complex<double>> zeros = abrth_roots(b);
        std::vector<std::complex<double>> poles = abrth_roots(a);

        if (!H_){   // if digital, invert poles and zeros
            for(int k = 0; k < poles.size(); k++){
                poles[k] = 1.0/poles[k];
                std::cout << poles[k] << "\n";
            }

            for(int k = 0; k < zeros.size(); k++){
                zeros[k] = 1.0/zeros[k];
                std::cout << zeros[k] << "\n";
            }
        }

        plt::figure();
        plt::set_aspect_equal();
        plt::title("PZ Plot");
        plt::grid(true);

        std::vector<double> real;
        std::vector<double> imag;

        // plot unit circle and axes
        double num_points = 200;
        for (int i = 0; i <= num_points; ++i) {
            double angle = 2.0 * M_PI * i / num_points;
            real.push_back(std::cos(angle));
            imag.push_back(std::sin(angle));
        }
        // plt::plot(real, imag, {{"color","#0072BD"}, {"linestyle", "--"}, {"linewidth", "0.8"}});
        plt::scatter(real, imag, 1.0, {{"color","#0072BD"},{"marker", "."}});
        real.clear();
        imag.clear();
    
        for (double i = 0; i <= num_points; ++i) {
            real.push_back((i/num_points-0.5)*6);
            imag.push_back(0);
        }
        plt::scatter(real, imag, 1.0, {{"color","#0072BD"},{"marker", "."}});
        real.clear();
        imag.clear();
    
        for (int i = 0; i <= num_points; ++i) {
            real.push_back(0);
            imag.push_back((i/num_points-0.5)*6);
        }
        plt::scatter(real, imag, 1.0, {{"color","#0072BD"},{"marker", "."}});
        real.clear();
        imag.clear();

        // plot zeros
        for (const auto& zero : zeros){
            real.push_back(zero.real());
            imag.push_back(zero.imag());
        }
        plt::scatter(real, imag, 40.0, {{"color", "blue"}, {"facecolor", "none"}}); 
        real.clear();
        imag.clear();

        // plot poles
        for (const auto& pole : poles){
            real.push_back(pole.real());
            imag.push_back(pole.imag());
        }
        plt::scatter(real, imag, 40.0, {{"color", "red"}, {"marker", "x"}}); 
        
        plt::xlabel("Real Part");
        plt::ylabel("Imaginary Part");

        plt::show();
    }

    void plot_ht(const std::vector<double>& b, const std::vector<double>& a, int N){
        std::vector<double> x_n(N, 0.0);
        x_n[0] = 1;

        plt::figure();
        plt::stem(filt(b, a, x_n));
        plt::title("Impulse Response h[n]");
        plt::xlabel("Sample");
        plt::ylabel("Magnitude");
        plt::grid(true);
    
        plt::show();
    }

    void plot_TF(const std::vector<double>& b, const std::vector<double>& a, TF_type H_, freq_type f, graph_type t, int N){
        std::vector<double> gain(N);
        std::vector<double> freqs(N);

        std::complex<double> num;
        std::complex<double> den;

        double FRQ = 100.0;

        for (double i = 0; i < N; i++){
            if (H_){
                freqs[i] = (f ? 1 : 2.0*M_PI)*(i/N)*FRQ;
            } else {
                freqs[i] = (f ? SAMPLE_RATE/2 : 1) * i/N;
            }

            num = 0;
            for (int n = 0; n < b.size(); n++){
                if (H_){
                    num += b[n] * std::pow(std::complex<double> (0.0, 2.0*M_PI*(i / N)*FRQ), n);
                } else {
                    num += std::polar(b[n], -n*(i/N)*M_PI);    // add b[n]*e^(+/-j*n*(i/N)*M_PI)
                }
            }

            den = 0;
            for (int d = 0; d < a.size(); d++){
                if (H_){
                    den += a[d] * std::pow(std::complex<double> (0.0, 2.0*M_PI*(i / N)*FRQ), d);
                } else {
                    den += std::polar(a[d], -d*(i/N)*M_PI);
                }
            }

            gain[i] = std::abs(num/den);
            if (!t){
                gain[i] = 20*std::log10(gain[i]);
            }
        }

        plt::figure();
        if (H_){
            plt::loglog(freqs, gain);
        } else {
            plt::plot(freqs, gain);
        }
        plt::title("H(z)");
        plt::xlabel((f ? "Frequency (Hz)" : (H_ ? "Frequency (Rad/s)" : "Normalized Frequency (x π rad/sample)")));
        plt::ylabel((t ? "Magnitude (abs)" : "Magnitude (dB)"));
        plt::grid(true);

        plt::show();
    }

    void plot_t(const std::vector<double>& x_n){
        int N = x_n.size();

        std::vector<double> n(N);

        for (double i = 0; i < N; i++){
            n[i] = i/SAMPLE_RATE;
        }

        plt::figure();
        plt::plot(n, x_n);
        plt::title("t-sig");
        plt::xlabel("Time (s)");
        plt::ylabel("Magnitude");
        plt::grid(true);

        plt::show();
    }


    std::vector<std::complex<double>> abrth_roots(const std::vector<double>& p, double tolerance, int maxIterations){
        std::vector<double> p_t = p;
        while (!p_t.empty() && p_t.back() == 0.0) { // pop trailing zeros
            p_t.pop_back();
        }
        int N = p_t.size() - 1;   // order of poly / number of roots
        if (!N){
            std::cerr << "Warning, zero size poly, no roots found\n\n";
            return std::vector<std::complex<double>> {};
        }
        std::vector<std::complex<double>> z(N);
        std::vector<std::complex<double>> w(N);

        // initialize equally spaced roots on unit circle
        for (int i = 0; i < N; i++){
            z[i] = std::polar(1.0, 2.0*M_PI*i/N);
        }

        for (int iter = 0; iter < maxIterations; iter++){

            // compute steps
            for (int k = 0; k < N; k++){
                
                // find p(x)
                std::complex<double> p_x = poly_eval(p_t, z[k]);

                // find p'(x)
                std::complex<double> p_x_prime = 0;
                for (int n = 1; n <= N; n++){
                    p_x_prime += n*p_t[n]*pow(z[k], n-1);
                }

                // find sum
                std::complex<double> sum = 0;
                for (int j = 0; j < N; j++){
                    if (j != k){
                        sum += 1.0/(z[k] - z[j]);
                    }
                }

                w[k] = (p_x/p_x_prime)/(1.0-(p_x/p_x_prime)*sum);
                // w[k] = p_x/p_x_prime;
            }

            // step roots
            bool converged = false;
            for (int k = 0; k < N; k++){
                z[k] -= w[k];
                if (std::abs(w[k]) < tolerance){
                    converged = true;
                }
            }

            if (converged) {
                std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
                return z;
            }
        }

        std::cout << "Aberth method did not converge within " << maxIterations << " iterations." << std::endl;
        return z;
    }

    void print_poly(const std::vector<double>& b, const char var, bool inv){
        for (int i = 0; i < b.size(); i++){
            if (b[i] != 0){
                if (i == 0) {
                    std::cout << b[i];
                } else {
                    std::cout << (b[i] >= 0 ? " + " : " - ");
                    if (abs(b[i] != 1)){
                        std::cout << abs(b[i]) << "*";
                    }
                    std::cout << var;
                    if (i != 1 || inv){
                        std::cout << "^" << (inv ? -i : i);
                    }
                }
            }
        }
        std::cout << std::endl;
    }

    void print_TF(const std::vector<double>& b, const std::vector<double>& a, TF_type H_){
        char var = 'z';
        bool inv = true;
        if (H_){    // analog settings
            var = 's';
            inv = false;
        }
        // print numenator
        print_poly(b, var, inv);
        std::cout << "------------------------------------------------------------------------\n";
        // print denominator
        print_poly(a, var, inv);
    }

    std::pair<std::vector<double>, std::vector<double>> butter(int N, double Omega_c, filter_type f){
        // Generate zeros and poles
        std::vector<std::complex<double>> zeros;
        if (f == LPF){ // no zeros
            zeros = {}; 
        } else if (f = HPF) {   // N repeated zeros
            zeros = std::vector<std::complex<double>> (N, 0.0); 
        }
        
        std::vector<std::complex<double>> poles;
        for (int i = 0; i < N; i++){     // find poles on LHS of S-plane
            poles.push_back(std::polar(Omega_c, ((2*i + 1)*M_PI)/(2.0*N) + M_PI_2));
        }

        auto results = zpk_to_ab(zeros, poles, std::pow(Omega_c, N));
        std::vector<double> num = results.first;
        std::vector<double> den = results.second;
        print_TF(num, den, ANA);

        return std::make_pair(num, den);
    }
    

    double butter_N(double WP, double Rp, double WS, double Rs, bool DB){
        Rs = (DB ? pow(10, Rs/10) : pow(1/Rs, 2));
        Rp = (DB ? pow(10, Rp/10) : pow(1/Rp, 2));
        double num = log10((Rp - 1)/(Rs - 1));
        double den = 2*log10(WP/WS);
        return num/den;
    }

    std::pair<std::vector<double>, std::vector<double>> zpk_to_ab(const std::vector<std::complex<double>>& z, const std::vector<std::complex<double>>& p, double k){
        // poles and zeros are in complex conjugate or purely real... only need to consider LHS of S-plane
        auto transformZ = [](const std::vector<std::complex<double>>& x) -> std::vector<double> {
            std::vector<double> y = {1};
            std::vector<double> term (3);
            int N = x.size();
            for (int i = 0; i < N/2.0; i++){     // find poles on LHS of S-plane
                // if ((N & 1) && (i == N-1)) {    // check if last pole on odd number of poles
                if (N/2.0 - i == 0.5) {    // check if last pole on odd number of poles
                    term[0] = -x[i].real();
                    term[1] = 1;
                    term[2] = 0;
                } else {
                    term[0] = pow(x[i].real(), 2) + pow(x[i].imag(), 2);
                    term[1] = -2*x[i].real();
                    term[2] = 1;
                }
                y = poly_multi(y, term);    // multiply new term into denominator
            }
            return y;
        };

        // N = z.size();
        // for (int i = 0; i < N; i++){     // find poles on LHS of S-plane
        //     if ((N & 1) && (i == N-1)) {    // check if last pole on odd number of poles
        //         term[0] = -z[i].real();
        //         term[1] = 1;
        //         term[2] = 0;
        //     } else {
        //         term[0] = pow(z[i].real(), 2) + pow(z[i].imag(), 2);
        //         term[1] = -2*z[i].real();
        //         term[2] = 1;
        //     }
        //     b = poly_multi(b, term);    // multiply new term into denominator
        // }
        return std::make_pair(poly_scale(transformZ(z), k), transformZ(p));
    }

}