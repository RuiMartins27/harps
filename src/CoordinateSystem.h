#ifndef COORDINATE_SYSTEM_H  // Include guard to prevent double inclusion
#define COORDINATE_SYSTEM_H

#include <iostream>
#include <string>
#include <complex>
#include <tuple>
#include <vector>
#include <unordered_set>
#include <array>

#include <petscksp.h>

#include "Constants.h"


using Complex = std::complex<double>;

class CoordinateSystem {
public:
    std::array<std::string, 3> fields_str = {"E_x", "E_y", "E_z"};
    std::vector<double> x, y, z;

    virtual ~CoordinateSystem() = default;
    
    virtual Mat createMaxwellEquationMatrix(Complex* f_grad_cond, Complex* complex_permittivity, double vacuum_wave_number, double waveguide_number, double yCenter, int size) = 0;

    virtual std::vector<Complex> calculateCondGradFunction(Complex* complex_conductivity, Complex* complex_permittivity, double angular_frequency) = 0;

    virtual void DirichletBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction , bool upper_boundary, Complex value, PetscInt low_rank, PetscInt high_rank) = 0;

    virtual void DirichletBoundaryConditionsPoint(Mat& matrix, Vec& b_array, Complex value, int idx, PetscInt low_rank, PetscInt high_rank) = 0;

    virtual void NeumannBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction , bool upper_boundary, Complex value, PetscInt low_rank, PetscInt high_rank) = 0;

    virtual void LocalNeumannCondition(Mat& matrix, Vec& b_array, Complex value, int idx_a, int idx_b, char direction, PetscInt low_rank, PetscInt high_rank) = 0;
    
    virtual void RobinBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction, bool upper_boundary, std::vector<Complex> wave_value, double vacuum_wave_number, double cutoff_wave_number, PetscInt low_rank, PetscInt high_rank) = 0;

    virtual void PeriodicBoundaryConditions(Mat& matrix, Vec& b_array, char direction, PetscInt low_rank, PetscInt high_rank) = 0;

    virtual void createPMLProfile(char direction, bool upper_boundary, int n_pml, double omega, int m, double sigma_0) = 0;

    virtual bool createNonUniformGrid(double refinementFactor_x, double refinementFactor_y, double refinementFactor_z, 
        double x_BL, double x_BR, double x_RL,double x_RR,
        double y_BL, double y_BR, double y_RL, double y_RR,
        double z_BL, double z_BR, double z_RL, double z_RR,
        std::vector<double> x_grid, std::vector<double> y_grid, std::vector<double> z_grid) = 0;

    virtual void AddMetalPointsInVolume(Mat& matrix, Vec& b_array, std::vector<std::tuple<int, int, int>> metalLocations, double vacuum_wave_number, PetscInt low_rank, PetscInt high_rank) = 0;
    
    virtual double* normalizeDens(double* dens, int num_points, double& total_abs_power) = 0;

    virtual void calculateFlux(Complex* fields, double* P_flux, double omega, int num_pts) = 0;

    virtual void calculatePoynting(double* P_flux, double* absorbedPowerDensity, double* P_poynting, int num_pts) = 0;
    
    virtual int point_index(int i, int j, int k=-1) = 0;

    virtual int getXindex(int idx) = 0;
    virtual int getYindex(int idx) = 0;
    virtual int getZindex(int idx) = 0;
    
    virtual double cellSize(const std::vector<double>& d, int idx, int n) = 0;
    virtual double getVolume(int i, int j, int k=-1) = 0;

    virtual double getX(int i) = 0;
    virtual double getY(int j) = 0;
    virtual double getZ(int k) = 0;

    virtual CoordinateSystem* get() = 0;


    double* computeAbsorbedPowerDens(const Complex* fields, double* real_conductivity, int num_points){
        double* P_abs_dens = new double[num_points];
        int idx;

        // P_abs_dens = 0.5 * Real{E dot J*} = 0.5 * Real{sigma} * |E|^2


        for (int i = 0; i < num_points; i++){
            idx = i*3; 
            P_abs_dens[i] = 0.5*real_conductivity[i]*(std::norm(fields[idx]) + std::norm(fields[idx+1]) + std::norm(fields[idx+2]));
        }

        return P_abs_dens;
    }

    void setRowToIdentity(Mat& matrix, int row) {

    };

    void setRowToZeroExcept(Mat& matrix, int row, int col, Complex value) {
      
    };

    void calculateFieldAmplitudes(const Complex* fields, double* E_amplitude, int num_points) {
        for(int point = 0; point < num_points; point++) {
            int idx = point * 3; 
            
            // Calculate |E| = sqrt(|Ex|^2 + |Ey|^2 + |Ez|^2)
            E_amplitude[point] = std::sqrt(std::norm(fields[idx]) + std::norm(fields[idx+1]) + std::norm(fields[idx+2]));
        }
    }
    void calculateFieldAmplitudes(const double* fields, double* E_amplitude, int num_points) {
        for(int point = 0; point < num_points; point++) {
            int idx = point * 3; 
            
            // Calculate |E| = sqrt(|Ex|^2 + |Ey|^2 + |Ez|^2)
            E_amplitude[point] = std::sqrt(fields[idx]*fields[idx] + fields[idx+1]*fields[idx+1] + fields[idx+2]*fields[idx+2]);
        }
    }

    void calculateFieldRealPart(const Complex* fields, double* E_real, int num_points){
        int idx;
        for(int point = 0; point < num_points; point++) {
            idx = point * 3; 
            
            // Calculate real{E} = sqrt((real{Ex})^2 + (real{Ey})^2 + |(real{Ez})^2)
            E_real[point] = std::sqrt(std::real(fields[idx])*std::real(fields[idx]) + std::real(fields[idx+1])*std::real(fields[idx+1]) + std::real(fields[idx+2])*std::real(fields[idx+2]));
        }
    }

    void calculateFieldImagPart(const Complex* fields, double* E_imag, int num_points){
        int idx;
        for(int point = 0; point < num_points; point++) {
            idx = point * 3; 
            
            // Calculate imag{E} = sqrt((imag{Ex})^2 + (imag{Ey})^2 + |(imag{Ez})^2)
            E_imag[point] = std::sqrt(std::imag(fields[idx])*std::imag(fields[idx]) + std::imag(fields[idx+1])*std::imag(fields[idx+1]) + std::imag(fields[idx+2])*std::imag(fields[idx+2]));
        }
    }

    void getFieldComponents(const Complex* fields, double* Ex_amp, double* Ey_amp, double* Ez_amp, int num_points) {
        for(int point = 0; point < num_points; point++) {
            int idx = point * 3;
            
            Ex_amp[point] = std::abs(fields[idx]);
            Ey_amp[point] = std::abs(fields[idx+1]);
            Ez_amp[point] = std::abs(fields[idx+2]);
        }
    }
    void getFieldComponents(const double* fields, double* Ex_amp, double* Ey_amp, double* Ez_amp, int num_points) {
        for(int point = 0; point < num_points; point++) {
            int idx = point * 3;
            
            Ex_amp[point] = fields[idx];
            Ey_amp[point] = fields[idx+1];
            Ez_amp[point] = fields[idx+2];
        }
    }

    virtual void calculatePlasmaFillingFactor(Complex* complex_conductivity, double y_center) = 0;
};




#endif // COORDINATE_SYSTEM_H
