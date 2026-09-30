#ifndef CARTESIAN_COORDINATES_H  // Include guard to prevent double inclusion
#define CARTESIAN_COORDINATES_H

#include "CoordinateSystem.h"

class CartesianCoordinateSystem : public CoordinateSystem {
private:
    int n_x, n_y, n_z;
    int n_points, n_variables;
    double length_x, length_y, length_z;      
    std::vector<double> dx, dy, dz;
    std::vector<Complex> s_profile_x, s_profile_y, s_profile_z;
    std::vector<Complex> dx_s, dy_s, dz_s;
public:
    CartesianCoordinateSystem(double num_x, double num_y, double num_z, double len_x, double len_y, double len_z) 
                            : n_x(num_x), n_y(num_y), n_z(num_z), length_x(len_x), length_y(len_y), length_z(len_z) {

        n_points = n_x * n_y * n_z;
        n_variables = 3 * n_points;

        for(int i=0; i<n_x; i++) x.push_back(i*length_x/(n_x-1)); 
        for(int j=0; j<n_y; j++) y.push_back(j*length_y/(n_y-1));
        for(int k=0; k<n_z; k++) z.push_back(k*length_z/(n_z-1));

        if(n_x == 1){
            x[0] = 0;
            dx.push_back(length_x);
        }
        if(n_y == 1){
            y[0] = 0;
            dy.push_back(length_y);
        }
        if(n_z == 1){
            z[0] = 0;
            dz.push_back(length_z);
        }

        for(int i=0; i<n_x-1; i++) dx.push_back((x[i+1]-x[i]));
        for(int j=0; j<n_y-1; j++) dy.push_back((y[j+1]-y[j]));
        for(int k=0; k<n_z-1; k++) dz.push_back((z[k+1]-z[k]));

        s_profile_x = std::vector<Complex>(n_x, Complex(1.0, 0.0));
        s_profile_y = std::vector<Complex>(n_y, Complex(1.0, 0.0));
        s_profile_z = std::vector<Complex>(n_z, Complex(1.0, 0.0));
    }
    
    Mat createMaxwellEquationMatrix(Complex* f_grad_cond, Complex* complex_permittivity, double vacuum_wave_number, double waveguide_number, double yCenter, int size);

    std::vector<Complex> calculateCondGradFunction(Complex* complex_conductivity, Complex* complex_permittivity, double angular_frequency);

    void DirichletBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction , bool upper_boundary, Complex value, PetscInt low_rank, PetscInt high_rank);

    void DirichletBoundaryConditionsPoint(Mat& matrix, Vec& b_array, Complex value, int idx, PetscInt low_rank, PetscInt high_rank);

    void NeumannBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction, bool upper_boundary, Complex value, PetscInt low_rank, PetscInt high_rank);

    void LocalNeumannCondition(Mat& matrix, Vec& b_array, Complex value, int idx_a, int idx_b, char direction, PetscInt low_rank, PetscInt high_rank);

    void RobinBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction, bool upper_boundary, std::vector<Complex> wave_value, double vacuum_wave_number, double cutoff_wave_number, PetscInt low_rank, PetscInt high_rank);

    void PeriodicBoundaryConditions(Mat& matrix, Vec& b_array, char direction, PetscInt low_rank, PetscInt high_rank);

    void createPMLProfile(char direction, bool upper_boundary, int n_pml, double omega, int m, double sigma_0)  override;

    bool createNonUniformGrid(double refinementFactor_x, double refinementFactor_y, double refinementFactor_z, 
                                double x_BL, double x_BR, double x_RL,double x_RR,
                                double y_BL, double y_BR, double y_RL, double y_RR,
                                double z_BL, double z_BR, double z_RL, double z_RR,
                                std::vector<double> x_grid, std::vector<double> y_grid, std::vector<double> z_grid);
    
    void AddMetalPointsInVolume(Mat& matrix, Vec& b_array, std::vector<std::tuple<int, int, int>> metalLocations, double y_reflector, PetscInt low_rank, PetscInt high_rank);

    void SetRowZeroManually(Mat& matrix, PetscInt row, PetscInt n_variables, int variables_per_point);

    double* normalizeDens(double* dens, int num_points, double& total_abs_power);

    void calculateFlux(Complex* fields, double* P_flux, double omega, int num_pts);

    void calculatePoynting(double* P_flux, double* absorbedPowerDensity, double* P_poynting, int num_pts);

    int point_index(int i, int j, int k) override { return 3*(i*n_y*n_z + j*n_z + k);}

    int getXindex(int idx) override { return (idx/3) / (n_y*n_z);}
    int getYindex(int idx) override { return (idx/3) / (n_z) % n_y;}
    int getZindex(int idx) override { return (idx/3) % n_z;}
    


    double cellSize(const std::vector<double>& d, int idx, int n) override {
        if (n == 1)         return d[0];
        if (idx == 0)       return 0.5*d[0];
        if (idx == n - 1)   return 0.5*d[n - 2];
        return 0.5*(d[idx] + d[idx - 1]);
    }

    double getVolume(int i, int j, int k) override{ return cellSize(dx, i, n_x)*cellSize(dy, j, n_y)*cellSize(dz, k, n_z);}

    double getX(int i) override { return x[i];}
    double getY(int j) override { return y[j];}
    double getZ(int k) override { return z[k];}

    void crossProduct(const Complex a_1, const Complex a_2, const Complex a_3, const Complex b_1, const Complex b_2, const Complex b_3,
                                                     Complex& r_1, Complex& r_2, Complex& r_3) {
        r_1 = a_2*b_3 - a_3*b_2; // x component
        r_2 = a_3*b_1 - a_1*b_3; // y component
        r_3 = a_1*b_2 - a_2*b_1; // z component
    }

    void calculatePlasmaFillingFactor(Complex* complex_conductivity, double y_center);

    virtual CoordinateSystem* get(){return this;}
};




#endif // CARTESIAN_COORDINATES_H
