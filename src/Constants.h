#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <complex>

using Complex = std::complex<double>;

extern std::string harps_dir;

namespace Constants {
    constexpr double EPSILON_0 = 8.854187817e-12;
    constexpr double MU_0 = 1.2566370614e-6;
    constexpr double C_LIGHT = 299792458.0;
    constexpr double CHARGE_E = 1.602176634e-19;
    constexpr double MASS_E = 9.1093837e-31;
    constexpr double K_BOLTZMANN = 1.380649e-23;
    constexpr Complex FREE{0.00123456789, -0.000123456789}; // Random but specific so it's easier to find in case of errors (signaling no excitation in this mode)
    constexpr Complex Z_0_I{0, MU_0 * C_LIGHT};
}

#endif // CONSTANTS_H
