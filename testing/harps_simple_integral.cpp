#include <complex>
#include <vector>
#include <cmath>
#include <memory>
#include <iostream>
#include <stdexcept>

using Complex = std::complex<double>;

namespace Constants {
    constexpr double EPSILON_0 = 8.854187817e-12;
    constexpr double MU_0 = 1.2566370614e-6;
    constexpr double C_LIGHT = 299792458.0;
    constexpr double CHARGE_E = 1.602176634e-19;
    constexpr double MASS_E = 9.1093837e-31;
    constexpr double K_BOLTZMANN = 1.380649e-23;
    constexpr Complex Z_0_I{0, MU_0 * C_LIGHT};
}


// This Code performs the WKB approximation to estimate the absorbed power in a plasma column

bool flag_uniform_field = 1;

int main(int argc, char* argv[]) {
    const Complex zero_C(0.0, 0.0);

    double R_in = 0.014;

    double total_abs_power = 0.0;
    
    // Get values from macro (niet)
    std::vector<double> ne, mu_real, mu_imag, radius;
    int n_pts = 200;

    for (int i = 0; i < n_pts; i++) {
        radius.push_back(i*R_in/(n_pts-1));
        ne.push_back(0.8e18 * exp(-radius[i]*radius[i] / (2*0.003*0.003)));  // in m^-3
        mu_real.push_back(0.8);
        mu_imag.push_back(2.5);
    }


    std::vector<double> real_permittivity(n_pts, 1.0);
    double angular_frequency = 2.0 * M_PI * 2.45e9;
    // double vacuum_wave_number = angular_frequency / Constants::C_LIGHT;

    double* real_conductivity = new double[n_pts];
    Complex* complex_conductivity = new Complex[n_pts];
    Complex* complex_permittivity = new Complex[n_pts];
    double* skin_depth = new double[n_pts];

    // Calculate plasma parameters
    for(int i = 0; i < n_pts; i++) {
        complex_conductivity[i] = Constants::CHARGE_E * ne[i] * (mu_real[i] + Complex(0,1)*mu_imag[i]);
        complex_permittivity[i] = real_permittivity[i] - Complex(0,1) * complex_conductivity[i] /  (angular_frequency * Constants::EPSILON_0);
        real_conductivity[i] = std::real(complex_conductivity[i]);

        // Components of k² = A + i*B 
        double k2_re = angular_frequency*angular_frequency*Constants::MU_0*Constants::EPSILON_0 - angular_frequency*Constants::MU_0*std::imag(complex_conductivity[i]);
        double k2_im = - angular_frequency*Constants::MU_0*std::real(complex_conductivity[i]);
                    
        // Full formula for Im(k)
        double k_imag = sqrt((sqrt(k2_re*k2_re + k2_im*k2_im) - k2_re) / 2.0);
        if(k2_im < 0) k_imag = -k_imag;

        skin_depth[i] = 1.0 / k_imag;  // m
        if(skin_depth[i] <= 0) skin_depth[i] = -skin_depth[i]; // choosing the positive root
    }

    double* E = new double[n_pts];
    E[n_pts-1] = 21000.0;  // arbitrary unit electric field at r=R_in
    for (size_t i = n_pts-1; i > 0; i--) {
        double dr = radius[i] - radius[i-1];
        if (flag_uniform_field) dr = 0;

        E[i-1] = E[i]*exp(-dr/(skin_depth[i-1]*0.25*R_in/radius[i-1]));  // exponential decay based on skin depth and coordinate transform
    }

    // Calculate absorbed power
    std::vector<double> P_abs(n_pts);
    for (int i = 0; i < n_pts; ++i) P_abs[i] = 0.5 * real_conductivity[i] * std::norm(E[i]);

    total_abs_power = 0.0;
    for (int i = 0; i < n_pts - 1; ++i) {
        double r1 = radius[i], r2 = radius[i + 1];
        double shell_area = M_PI * (r2 * r2 - r1 * r1);
        total_abs_power += 0.5 * (P_abs[i] + P_abs[i + 1]) * shell_area;
    }

    // Export results to macro
    for (int i = 0; i < n_pts; ++i)  std::cout << radius[i] << ", " << sqrt(std::norm(E[i])) <<  ",  " << P_abs[i] <<  std::endl;

    delete[] real_conductivity; delete[] complex_conductivity; delete[] complex_permittivity; delete[] skin_depth; delete[] E;
}