#ifndef CONFIG_PARSER_H  // Include guard to prevent double inclusion
#define CONFIG_PARSER_H

#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <regex>
#include <vector>
#include <tuple>

#include "Constants.h"


using Complex = std::complex<double>;

struct HARPSConfig {
    // Grid parameters
    int n_x, n_y, n_z = 1;
    double lengthX, lengthY, lengthZ;
    double refinementFactor_x = -1; double refinementFactor_y = -1; double refinementFactor_z = -1;
    double x_BL, x_BR, x_RL, x_RR; // buffer left position, buffer right position, refinement left position, refinement right position 
    double y_BL, y_BR, y_RL, y_RR;
    double z_BL, z_BR, z_RL, z_RR;
    std::vector<double> x_grid;
    std::vector<double> y_grid;
    std::vector<double> z_grid;
    
    // Physical parameters
    double frequency;
    double waveguide_number = M_PI/0.08636; // default waveguide number corresponds to TE10 on WR340
    double yCenter = 0.056;
    double realPermittivity;
    std::vector<std::tuple<int, int, int>> realPermittivityLocations;
    std::vector<double> realPermittivityValues;
    double electronDensity;
    std::vector<std::tuple<int, int, int>> electronDensityLocations;
    std::vector<double> electronDensityValues;
    double realMobility;
    std::vector<std::tuple<int, int, int>> realMobilityLocations;
    std::vector<double> realMobilityValues;
    double imagMobility;
    std::vector<std::tuple<int, int, int>> imagMobilityLocations;
    std::vector<double> imagMobilityValues;
    
    // Solver parameters
    std::string solverMethod;
    std::string PC = "jacobi";
    int scalar = -1;
    double tolerance;
    int maxIterations;
    std::string initial_guess_file = "";
    bool flag_cylindrical_plasma = 1;
    
    // Boundary conditions
    std::string xLowerBC, xUpperBC;
    std::string yLowerBC, yUpperBC;
    std::string zLowerBC, zUpperBC;

    std::vector<Complex> excitationValue;

    Complex DirichletBx_x = Complex(0,0);
    Complex DirichletBx_y = Complex(1,0);
    Complex DirichletBx_z = Complex(1,0);

    Complex DirichletBy_x = Complex(1,0);
    Complex DirichletBy_y = Complex(0,0);
    Complex DirichletBy_z = Complex(1,0);

    Complex DirichletBz_x = Complex(1,0);
    Complex DirichletBz_y = Complex(1,0);
    Complex DirichletBz_z = Complex(0,0);

    std::vector<std::vector<Complex>> injectionValues;
    
    std::vector<std::tuple<int, int, int>> sourceLocations;
    std::vector<std::vector<Complex>> sourceValues;
    std::vector<std::tuple<int, int, int>> excitationLocations;
    std::vector<std::vector<Complex>> excitationValues;
    std::vector<std::tuple<int, int, int>> metalLocations;

    double y_reflector = 1e6; // default value for y_reflector is set to a large number, effectively disabling it unless specified on input file
    
    // PML parameters
    bool enableXLowerPML = false, enableXUpperPML = false, enableYLowerPML = false, enableYUpperPML = false;
    bool enableZLowerPML = false, enableZUpperPML = false;
    int xLowerLayers = 0, xUpperLayers = 0, yLowerLayers = 0, yUpperLayers = 0;
    int zLowerLayers = 0, zUpperLayers = 0;
    int orderPML;
    double sigma_0 = 0.8;
    
    // Output parameters
    bool printConfig = false;
    bool printGrid = false;
    bool printProgress = false;
    bool printMatrix = false;
    bool printSolution = false;
    bool storeResults = false;
    std::string outputDirectory;
    std::vector<std::string> outputFields;
};

class ConfigParser {
public:
    ConfigParser() = default;

    // Reads an input file and returns the parsed HARPSConfig.
    HARPSConfig parseFile(const std::string& filename);

private:
    std::string trim(const std::string& str);
    std::vector<std::string> split(const std::string& str, char delimiter);
    std::vector<std::tuple<int, int, int>> parseLocationArray(const std::string& str, const HARPSConfig& config);
    std::vector<std::vector<Complex>> parseComplexValueArray(const std::string& str, const HARPSConfig& config);
    std::vector<Complex> parseComplexValueRow(const std::string& str);
    std::vector<double> parseDoubleValueArray(const std::string& str, const HARPSConfig& config);
    int parseExpression(const std::string& expr, const HARPSConfig& config);

    // Stream-based file parsers
    std::vector<std::vector<Complex>> parseComplexValueArrayFromFile(const std::string& filename);
    std::vector<std::tuple<int, int, int>> parseLocationArrayFromFile(const std::string& filename, const HARPSConfig& config);
    std::vector<double> parseDoubleValueArrayFromFile(const std::string& filename);
    
    std::tuple<int, int, int> parseLocationTuple(const std::string& pointStr, const HARPSConfig& config);
    int string_to_int(std::string str, std::string key="internal conversions from input file", char cond='a');
    double string_to_double(std::string str, std::string key="internal conversions from input file", char cond='a');
    bool parseBoolean(const std::string& value, std::string key);
};

// Non-member operator <<

std::ostream& operator<<(std::ostream& os, const std::tuple<int, int, int>& t);

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    for (const auto& v : vec) {
        os << v << " ";
    }
    os << "]";
    return os;
}


inline std::ostream& operator<<(std::ostream& os, const std::vector<Complex>& vec) {
    os << "[ ";
    for (const auto& v : vec) {
        if (std::abs(v.real() - Constants::FREE.real()) < 1e-10 && std::abs(v.imag() - Constants::FREE.imag()) < 1e-10) {
            os << "F ";
        } else {
            os << v << " ";
        }
    }
    os << "]";
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<Complex>>& vec) {
    os << "[ ";
    for (const auto& inner_vec : vec) {
        os << "[ ";
        for (const auto& v : inner_vec) {
            if (std::abs(v.real() - Constants::FREE.real()) < 1e-10 && std::abs(v.imag() - Constants::FREE.imag()) < 1e-10) {
                os << "F ";
            } else {
                os << v << " ";
            }
        }
        os << "] ";
    }
    os << "]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const HARPSConfig& config);


#endif // CONFIG_PARSER_H
