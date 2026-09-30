#include "ConfigParser.h"
#include "CartesianCoordinates.h"
#include "MatrixSolver.h"
#include "Constants.h"
#include "Utils.h"
#include <cstring>

#include <libmuscle/libmuscle.hpp>
#include <ymmsl/ymmsl.hpp>


std::string harps_dir = "../../../../harps/";    // Path to the harps directory if I'm running from harps or 1D_radial if it's next to harps


using libmuscle::Data;
using libmuscle::DataConstRef;
using libmuscle::Instance;
using libmuscle::Message;
using libmuscle::PortsDescription;

using ymmsl::Operator;


bool receive_data(Instance & client, const std::string & port, bool & flag_out, std::vector<double> & vec_out, double & dt) {
    Message msg = client.receive(port);
    auto data = msg.data();
    dt = msg.timestamp();

    if (data.is_a<bool>()) {
        flag_out = data.as<bool>();

        std::cout << "Boolean received on port " << port << ": " << std::boolalpha << flag_out << std::endl;
        return false;  // Flag
    }

    else{
        vec_out.clear();
        auto data_ptr = data.elements<double>();
        std::vector<double> U(data_ptr, data_ptr + data.size());
        vec_out = std::move(U);

        std::cout << "List received on port " << port << ": " << vec_out[0] << ", " << vec_out[1] << ", ..." << std::endl;
        return true;  // Vector
    }

    throw std::runtime_error("Unsupported message type received");
}

double angle(double x, double y){ // [0, 2pi];  being (0,-1) direction since it's incoming wave region
    if(abs(x) < 1e-10){
        if(y > 0) return M_PI;
        else return 0.0;
    }
    
    if(x>0) return atan(y/x) + M_PI/2;
    else    return atan(y/x) + 3*M_PI/2;
}


void harps_main(int argc, char* argv[]){
    using Complex = std::complex<double>;
    const Complex zero_C(0.0, 0.0);

    int rank, size;
    const char* const* argv_const = const_cast<const char* const*>(argv);
    PetscInitialize(&argc, &argv, NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    const bool symmetric = true; // Symmetric harps over center x=C_x plane
    const bool symmetric_harps = true;
    const int N_macro = 8;
    const double T_e = 10000.0; // Assumed electron temperature in Kelvin

    double angle_macro_2;
    if(symmetric) angle_macro_2 = 0.5*M_PI/(N_macro-1);
    else angle_macro_2 = M_PI/N_macro;

    std::unique_ptr<Instance> instance;

    if (rank == 0) {
        std::vector<std::string> beginHarps;
        std::vector<std::string> endHarps;

        for (int i = 0; i < N_macro; ++i) {
            beginHarps.push_back("begin_harps_" + std::to_string(i));
            endHarps.push_back("end_harps_" + std::to_string(i));
        }

        PortsDescription ports{ {Operator::F_INIT, beginHarps}, {Operator::O_F, endHarps} };

        instance = std::make_unique<Instance>(argc, argv_const, ports);
    }
    
    std::vector<double> r_bins, r_centers, radial_power, radial_volume, U;
    std::vector<std::vector<double>> r_P_ne_out(N_macro);
    double x, y, r, theta, volume, pabs;
    bool do_harps, flag_sync;

    double center_x = 0.04318;
    if(symmetric_harps) center_x = 0.0;
    double center_y = 0.046;
    double R_in = 0.014;

    bool instance_cond = true;


    while (instance_cond) {
        do_harps = false;
        flag_sync = false;
        double dt = 0.0;
        std::array<std::vector<double>, N_macro> data_in;

        if(rank == 0){
            instance_cond = instance->reuse_instance();

            for(int i = 0; i < N_macro; ++i) {
                bool flag_out = false;
                double dt_out;

                do_harps = receive_data(*instance, "begin_harps_" + std::to_string(i), flag_out, data_in[i], dt_out);
                flag_sync = flag_out || flag_sync;
                if (dt_out > dt) dt = dt_out;
            }
        }
        if (dt < 1e-10) dt = 0.5e-6; // Default time step if not set

        std::cout << dt << " s time step received" << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&do_harps, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&instance_cond, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

       double total_abs_power = 0.0;

        // real harps things
        if(do_harps){
            try {
                ConfigParser parser;
                HARPSConfig config;
                
                // Start measuring time
                PetscLogDouble startTime, endTime;
                PetscTime(&startTime);
                double previous_time = startTime;

                config = parser.parseFile(harps_dir  + argv[1]);

                if(config.printConfig && rank == 0) std::cout << config << std::endl;

                // debugging prints
                bool print_ranks = 0;
                bool print_times = 0;

                int num_pts, num_variables, temp_idx;
                double angular_frequency, vacuum_wave_number;

                num_pts = (config.n_x * config.n_y * config.n_z);
                num_variables = num_pts * 3;

                if(num_pts%size != 0){
                    PetscFinalize();
                    throw std::runtime_error("Number of threads (" + std::to_string(size) +") needs to evenly divide the number of points (" + std::to_string(num_pts) +"): " + std::to_string(1.*num_pts/size));
                }
                
                angular_frequency = 2.0 * M_PI * config.frequency;
                vacuum_wave_number = angular_frequency / Constants::C_LIGHT;

                // Initialize arrays
                double* electron_density = new double[num_pts];
                double* real_permittivity = new double[num_pts];
                double* real_conductivity = new double[num_pts];
                double* real_mobility = new double[num_pts];
                double* imag_mobility = new double[num_pts];
                Complex* complex_conductivity = new Complex[num_pts];
                Complex* complex_permittivity = new Complex[num_pts];




                std::unique_ptr<CoordinateSystem> Grid;
                Grid = std::make_unique<CartesianCoordinateSystem>(config.n_x, config.n_y, config.n_z, config.lengthX, config.lengthY, config.lengthZ);
                bool non_uniform_grid = Grid->createNonUniformGrid(config.refinementFactor_x, config.refinementFactor_y, config.refinementFactor_z, 
                                    config.x_BL, config.x_BR, config.x_RL, config.x_RR, config.y_BL, config.y_BR, config.y_RL, config.y_RR,
                                    config.z_BL, config.z_BR, config.z_RL, config.z_RR, config.x_grid, config.y_grid, config.z_grid);

                if(config.printGrid && rank == 0 && non_uniform_grid) printGrid(config.outputDirectory + "Grid.txt", config.n_x, config.n_y, config.n_z,
                                                                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);


                // Initialize default values
                std::fill_n(electron_density, num_pts, config.electronDensity);
                std::fill_n(real_permittivity, num_pts, config.realPermittivity);
                std::fill_n(real_mobility, num_pts, config.realMobility);
                std::fill_n(imag_mobility, num_pts, config.imagMobility);

                // Apply input position specific values of epsilon_r and n_e

                // Fill in the n_e and nu_eh arrays with input from macro (might need to interpolate)!!!!!!!!!!!!!!
                // Assume U = {rad_1, ne_1, mu_real_1, mu_imag_1, rad_2, ne_2, mu_real_2, mu_imag_2, ...}

                if(rank ==0){
                    std::array<std::vector<double>, N_macro> ne, mu_real, mu_imag, radius;

                    for (int i = 0; i < N_macro; ++i) {
                        for (std::size_t j = 0; j < data_in[i].size(); j += 4) {
                            radius[i].push_back(data_in[i][j]);
                            ne[i].push_back(data_in[i][j+1]);
                            mu_real[i].push_back(data_in[i][j+2]);
                            mu_imag[i].push_back(data_in[i][j+3]);
                        }

                        radius[i].pop_back(); radius[i].pop_back(); ne[i].pop_back(); ne[i].pop_back(); mu_real[i].pop_back(); mu_real[i].pop_back(); mu_imag[i].pop_back(); mu_imag[i].pop_back(); // Remove last two points to avoid extrapolation
                    }

                    std::vector<MonotoneCubicInterpolation<double>> spline_sets_ne;
                    std::vector<MonotoneCubicInterpolation<double>> spline_sets_mu_real;
                    std::vector<MonotoneCubicInterpolation<double>> spline_sets_mu_imag;
                    for (long unsigned int i = 0; i < radius[0].size(); ++i) {
                        double step;
                        std::vector<double> x, y_ne, y_mu_real, y_mu_imag;

                        if(symmetric){
                            step = M_PI/(N_macro-1);

                            x.push_back(-2*step); y_ne.push_back(ne[2][i]); y_mu_real.push_back(mu_real[2][i]); y_mu_imag.push_back(mu_imag[2][i]);
                            x.push_back(-step);   y_ne.push_back(ne[1][i]); y_mu_real.push_back(mu_real[1][i]); y_mu_imag.push_back(mu_imag[1][i]);

                            // Points from macros
                            for(int j = 0; j < N_macro; ++j){
                                x.push_back(j*step);
                                y_ne.push_back(ne[j][i]);
                                y_mu_real.push_back(mu_real[j][i]);
                                y_mu_imag.push_back(mu_imag[j][i]);
                            }
                            // Points from macros reflected
                            for(int j = N_macro-2; j > 0; --j){
                                x.push_back((2*(N_macro-1) - j)*step);
                                y_ne.push_back(ne[j][i]);
                                y_mu_real.push_back(mu_real[j][i]);
                                y_mu_imag.push_back(mu_imag[j][i]);
                            }

                            x.push_back(2*M_PI);        y_ne.push_back(ne[0][i]); y_mu_real.push_back(mu_real[0][i]); y_mu_imag.push_back(mu_imag[0][i]);
                            x.push_back(2*M_PI + step); y_ne.push_back(ne[1][i]); y_mu_real.push_back(mu_real[1][i]); y_mu_imag.push_back(mu_imag[1][i]);

                            spline_sets_ne.emplace_back(x.data(), y_ne.data(), x.size());
                            spline_sets_mu_real.emplace_back(x.data(), y_mu_real.data(), x.size());
                            spline_sets_mu_imag.emplace_back(x.data(), y_mu_imag.data(), x.size());
                        } else{
                            step = 2*M_PI/N_macro;

                            x.push_back(-2*step); y_ne.push_back(ne[N_macro-2][i]); y_mu_real.push_back(mu_real[N_macro-2][i]); y_mu_imag.push_back(mu_imag[N_macro-2][i]); // periodic points
                            x.push_back(-step);   y_ne.push_back(ne[N_macro-1][i]); y_mu_real.push_back(mu_real[N_macro-1][i]); y_mu_imag.push_back(mu_imag[N_macro-1][i]);

                            // Points from macros
                            for(int j = 0; j < N_macro; ++j){
                                x.push_back(j*step);
                                y_ne.push_back(ne[j][i]);
                                y_mu_real.push_back(mu_real[j][i]);
                                y_mu_imag.push_back(mu_imag[j][i]);
                            }
                            x.push_back(2*M_PI);        y_ne.push_back(ne[0][i]); y_mu_real.push_back(mu_real[0][i]); y_mu_imag.push_back(mu_imag[0][i]);
                            x.push_back(2*M_PI + step); y_ne.push_back(ne[1][i]); y_mu_real.push_back(mu_real[1][i]); y_mu_imag.push_back(mu_imag[1][i]);

                            spline_sets_ne.emplace_back(x.data(), y_ne.data(), x.size());
                            spline_sets_mu_real.emplace_back(x.data(), y_mu_real.data(), x.size());
                            spline_sets_mu_imag.emplace_back(x.data(), y_mu_imag.data(), x.size());
                        }
                    }

                    for (int i = 0; i < num_pts; ++i) {
                        x = Grid->getX(Grid->getXindex(i*3));
                        y = Grid->getY(Grid->getYindex(i*3));
                        r = std::sqrt((x - center_x)*(x - center_x) + (y - center_y)*(y - center_y));
                        theta = angle(x - center_x, y - center_y);

                        if (r > radius[0][radius[0].size()-1]) continue;
                        
                        std::pair<int,int> r_idx = find_indices(radius[0], r);

                        electron_density[i] = linear_interpolation_binary(radius[0][r_idx.first], spline_sets_ne[r_idx.first](theta), radius[0][r_idx.second], spline_sets_ne[r_idx.second](theta), r);
                        real_mobility[i] = linear_interpolation_binary(radius[0][r_idx.first], spline_sets_mu_real[r_idx.first](theta), radius[0][r_idx.second], spline_sets_mu_real[r_idx.second](theta), r);
                        imag_mobility[i] = linear_interpolation_binary(radius[0][r_idx.first], spline_sets_mu_imag[r_idx.first](theta), radius[0][r_idx.second], spline_sets_mu_imag[r_idx.second](theta), r);
                    }


                    // Apply Diffusion on ne
                    double* electron_density_old = new double[num_pts];
                    for (int i = 0; i < num_pts; ++i)   electron_density_old[i] = electron_density[i]; 

                    double dr = 0.00025;
                    int num_bins = static_cast<int>(std::ceil(R_in / dr));

                    std::vector<std::vector<int>> bin_point_indices(num_bins);
                    std::vector<std::vector<double>> bin_point_x(num_bins), bin_point_y(num_bins);
                    std::vector<std::vector<double>> bin_point_theta(num_bins);
                    std::vector<std::vector<double>> bin_point_volume(num_bins);

                    // --- Bin points by radius ---
                    for (int i = 0; i < config.n_x; ++i) {
                        for (int j = 0; j < config.n_y; ++j) {
                            int idx = Grid->point_index(i, j) / 3;
                            double x = Grid->getX(i);
                            double y = Grid->getY(j);
                            double dx = x - center_x;
                            double dy = y - center_y;
                            double r = std::sqrt(dx * dx + dy * dy);
                            double theta = angle(dx, dy);
                            double volume = Grid->getVolume(i, j, 0);

                            int bin = static_cast<int>(r / dr);
                            if (bin < num_bins){
                                bin_point_indices[bin].push_back(idx);
                                bin_point_x[bin].push_back(x);
                                bin_point_y[bin].push_back(y);
                                bin_point_theta[bin].push_back(theta);
                                bin_point_volume[bin].push_back(volume);
                            }
                        }
                    }
                    
                    // Convolution with a Gaussian kernel
                    for (int bin = 0; bin < num_bins; ++bin) {
                        r = (bin+0.5)*dr;
                        int bin_size = bin_point_indices[bin].size();

                        for (int j = 0; j < bin_size; ++j) {
                            double Tg = 4000;
                            double average_mu = 2.2e-4;
                            int idx = bin_point_indices[bin][j];
                            double theta = bin_point_theta[bin][j];

                            double D_a = (Tg+T_e)*Constants::K_BOLTZMANN*average_mu/(r*r*Constants::CHARGE_E); // rad^2/s
                            double sigma_D = std::sqrt(2*D_a*dt); // rad

                            double inv_sigma = 1.0/sigma_D;
                            double inv_sigma2 = inv_sigma*inv_sigma;
                            double norm_factor = inv_sigma/std::sqrt(2 * M_PI);

                            double sum = 0;
                            double norm = 0;
                            for (int m = 0; m < bin_size; ++m) {
                                int idx_m = bin_point_indices[bin][m];
                                double dtheta = bin_point_theta[bin][m] - theta;
                                
                                if (dtheta > M_PI) dtheta -= 2*M_PI;
                                else if (dtheta < -M_PI) dtheta += 2*M_PI;

                                if(abs(dtheta) > 3*sigma_D) continue;

                                double G = std::exp(-dtheta*dtheta*0.5*inv_sigma2)*norm_factor;

                                sum += electron_density_old[idx_m]*G*bin_point_volume[bin][m];
                                norm += G*bin_point_volume[bin][m];
                            }
                            if(norm > 1e-20) electron_density[idx] = sum / norm;
                        }
                    }
                    delete[] electron_density_old;
                }

                // Broadcast the values received in rank 0 to all others
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Bcast(electron_density, num_pts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(real_mobility, num_pts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast(imag_mobility, num_pts, MPI_DOUBLE, 0, MPI_COMM_WORLD);


                for (size_t i = 0; i < config.electronDensityLocations.size(); ++i) {
                    const auto& loc = config.electronDensityLocations[i];
                    double value = config.electronDensityValues[i];
                    
                    size_t index = std::get<0>(loc)*config.n_z*config.n_y + std::get<1>(loc)*config.n_z + std::get<2>(loc);
                    electron_density[index] = value;
                }
                for (size_t i = 0; i < config.realPermittivityLocations.size(); ++i) {
                    const auto& loc = config.realPermittivityLocations[i];
                    double value = config.realPermittivityValues[i];
                    
                    size_t index = std::get<0>(loc)*config.n_z*config.n_y + std::get<1>(loc)*config.n_z + std::get<2>(loc);
                    real_permittivity[index] = value;
                }
                for (size_t i = 0; i < config.realMobilityLocations.size(); ++i) {
                    const auto& loc = config.realMobilityLocations[i];
                    double value = config.realMobilityValues[i];
                    
                    size_t index = std::get<0>(loc)*config.n_z*config.n_y + std::get<1>(loc)*config.n_z + std::get<2>(loc);
                    real_mobility[index] = value;
                }
                for (size_t i = 0; i < config.imagMobilityLocations.size(); ++i) {
                    const auto& loc = config.imagMobilityLocations[i];
                    double value = config.imagMobilityValues[i];
                    
                    size_t index = std::get<0>(loc)*config.n_z*config.n_y + std::get<1>(loc)*config.n_z + std::get<2>(loc);
                    imag_mobility[index] = value;
                }

                // Create output directory if it doesn't exist
                config.outputDirectory = harps_dir + config.outputDirectory;
                if (!config.outputDirectory.empty() && config.outputDirectory.back() != '/') config.outputDirectory += '/';
                
                if(print_times){
                    PetscTime(&endTime);
                    PetscPrintf(PETSC_COMM_WORLD, "Time taken on config: %.3f seconds\n", endTime - previous_time);
                    previous_time = endTime;
                }

                if(config.printProgress && rank == 0) std::cout << "Creating System Matrix (" << num_variables << "x" << num_variables << ")" << std::endl;

                // Apply PML if enabled (needs to be called before creating matrix)
                if (config.enableXLowerPML)  Grid->createPMLProfile('x', false, config.xLowerLayers, angular_frequency, config.orderPML, config.sigma_0);
                if (config.enableXUpperPML)  Grid->createPMLProfile('x', true, config.xUpperLayers, angular_frequency, config.orderPML, config.sigma_0);
                if (config.enableYLowerPML)  Grid->createPMLProfile('y', false, config.yLowerLayers, angular_frequency, config.orderPML, config.sigma_0);
                if (config.enableYUpperPML)  Grid->createPMLProfile('y', true, config.yUpperLayers, angular_frequency, config.orderPML, config.sigma_0);
                if (config.enableZLowerPML)  Grid->createPMLProfile('z', false, config.zLowerLayers, angular_frequency, config.orderPML, config.sigma_0);
                if (config.enableZUpperPML)  Grid->createPMLProfile('z', true, config.zUpperLayers, angular_frequency, config.orderPML, config.sigma_0);

                // Calculate plasma conductivity and permitivity
                for (int i = 0; i < num_pts; ++i) {
                    Complex mobility(real_mobility[i], imag_mobility[i]);
                    complex_conductivity[i] = Constants::CHARGE_E*electron_density[i]*mobility;
                }
                if (config.flag_cylindrical_plasma) Grid->calculatePlasmaFillingFactor(complex_conductivity, config.yCenter);     // 2D YZ Plasma Filling in X

                for (int i = 0; i < num_pts; ++i){
                    real_conductivity[i] = std::real(complex_conductivity[i]);  // Used for p_abs = 0.5*cond_real*|E|^2
                    complex_permittivity[i] = real_permittivity[i] - Complex(0.0, 1.0)*complex_conductivity[i]/(angular_frequency*Constants::EPSILON_0);
                }


                Complex* f_grad_cond = Grid->calculateCondGradFunction(complex_conductivity, complex_permittivity, angular_frequency);
                
                // Create system matrix
                Mat systemMaxwell = Grid->createMaxwellEquationMatrix(f_grad_cond, complex_permittivity, vacuum_wave_number, config.waveguide_number, size);
                
                Vec b_vector;
                VecCreate(PETSC_COMM_WORLD, &b_vector);
                VecSetSizes(b_vector, PETSC_DECIDE, num_variables);
                VecSetFromOptions(b_vector);
                VecSetBlockSize(b_vector, 3);

                // Set right-hand side (forcing function)
                VecSet(b_vector, 0.0);

                PetscInt rstart, rend;
                MatGetOwnershipRange(systemMaxwell, &rstart, &rend);
                PetscInt low_rank, high_rank;
                VecGetOwnershipRange(b_vector, &low_rank, &high_rank);
                
                if(print_ranks){
                    std::cout << "Vector range: " << rank << " owns rows " << low_rank << " to " << high_rank << std::endl;
                    std::cout << "Matrix range: " << rank << " owns rows " << rstart   << " to " << rend << std::endl;
                }

                MatAssemblyBegin(systemMaxwell, MAT_FLUSH_ASSEMBLY);
                MatAssemblyEnd(systemMaxwell, MAT_FLUSH_ASSEMBLY);

                if(print_times){
                    PetscTime(&endTime);
                    PetscPrintf(PETSC_COMM_WORLD, "Time taken on creating Matrix: %.3f seconds\n", endTime - previous_time);
                    previous_time = endTime;
                }

                // Apply boundary conditions
                if(config.printProgress  && rank == 0) std::cout << "Applying Boundary Conditions" << std::endl;

                auto applyXBC = [&](const std::string& bcType, bool upper) {
                    if (bcType == "PerfectConductor") {
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "E_y", 'x', upper, zero_C, low_rank, high_rank);
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "E_z", 'x', upper, zero_C, low_rank, high_rank);
                        Grid->NeumannBoundaryConditions(systemMaxwell, b_vector, "E_x", 'x', upper, zero_C, low_rank, high_rank);
                    } else if (bcType == "Homogeneous") {
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "all", 'x', upper, zero_C, low_rank, high_rank);
                    } else if (bcType == "Neumann") {
                        Grid->NeumannBoundaryConditions(systemMaxwell, b_vector, "all", 'x', upper, zero_C, low_rank, high_rank);
                    }  else if (bcType == "Robin") {
                        for(int j = 0; j<3; j++) Grid->RobinBoundaryConditions(systemMaxwell, b_vector, Grid->fields_str[j], 'x', upper, config.injectionValues[j], vacuum_wave_number, M_PI/0.08636, low_rank, high_rank);
                    } else if (bcType == "Excitation") {
                        for(int j = 0; j<3; j++){
                            if (std::abs(config.excitationValue[j].real() - Constants::FREE.real()) > 1e-10 && std::abs(config.excitationValue[j].imag() - Constants::FREE.imag()) > 1e-10)
                                Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, Grid->fields_str[j], 'x', upper, config.excitationValue[j], low_rank, high_rank);
                        }
                    }
                };

                auto applyYBC = [&](const std::string& bcType, bool upper) {
                    if (bcType == "Robin") {
                        for(int j = 0; j<3; j++) Grid->RobinBoundaryConditions(systemMaxwell, b_vector, Grid->fields_str[j], 'y', upper, config.injectionValues[j], vacuum_wave_number, M_PI/0.08636, low_rank, high_rank);
                    }else if (bcType == "PerfectConductor") {
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "E_x", 'y', upper, zero_C, low_rank, high_rank);
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "E_z", 'y', upper, zero_C, low_rank, high_rank);
                        Grid->NeumannBoundaryConditions(systemMaxwell, b_vector, "E_y", 'y', upper, zero_C, low_rank, high_rank);
                    } else if (bcType == "Homogeneous") {
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "all", 'y', upper, zero_C, low_rank, high_rank);
                    } else if (bcType == "Neumann") {
                        Grid->NeumannBoundaryConditions(systemMaxwell, b_vector, "all", 'y', upper, zero_C, low_rank, high_rank);
                    } else if (bcType == "Excitation") {
                        for(int j = 0; j<3; j++){
                            if (std::abs(config.excitationValue[j].real() - Constants::FREE.real()) > 1e-10 && std::abs(config.excitationValue[j].imag() - Constants::FREE.imag()) > 1e-10)
                                Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, Grid->fields_str[j], 'y', upper, config.excitationValue[j], low_rank, high_rank);
                        }
                    }
                };

                auto applyZBC = [&](const std::string& bcType, bool upper) {
                    if (bcType == "PerfectConductor") {
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "E_x", 'z', upper, zero_C, low_rank, high_rank);
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "E_y", 'z', upper, zero_C, low_rank, high_rank);
                        Grid->NeumannBoundaryConditions(systemMaxwell, b_vector, "E_z", 'z', upper, zero_C, low_rank, high_rank);
                    } else if (bcType == "Homogeneous") {
                        Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, "all", 'z', upper, zero_C, low_rank, high_rank);
                    }  else if (bcType == "Neumann") {
                        Grid->NeumannBoundaryConditions(systemMaxwell, b_vector, "all", 'z', upper, zero_C, low_rank, high_rank);
                    }  else if (bcType == "Robin") {
                        for(int j = 0; j<3; j++) Grid->RobinBoundaryConditions(systemMaxwell, b_vector, Grid->fields_str[j], 'z', upper, config.injectionValues[j], vacuum_wave_number, M_PI/0.08636, low_rank, high_rank);
                    } else if (bcType == "Excitation") {
                        for(int j = 0; j<3; j++){
                            if (std::abs(config.excitationValue[j].real() - Constants::FREE.real()) > 1e-10 && std::abs(config.excitationValue[j].imag() - Constants::FREE.imag()) > 1e-10)
                                Grid->DirichletBoundaryConditions(systemMaxwell, b_vector, Grid->fields_str[j], 'z', upper, config.excitationValue[j], low_rank, high_rank);
                        }
                    }
                };
                
                applyXBC(config.xLowerBC, false); applyXBC(config.xUpperBC, true);
                applyYBC(config.yLowerBC, false); applyYBC(config.yUpperBC, true);
                applyZBC(config.zLowerBC, false); applyZBC(config.zUpperBC, true);

                if(config.xUpperBC == "Periodic") Grid->PeriodicBoundaryConditions(systemMaxwell, b_vector,'x', low_rank, high_rank);
                if(config.yUpperBC == "Periodic") Grid->PeriodicBoundaryConditions(systemMaxwell, b_vector,'y', low_rank, high_rank);
                if(config.zUpperBC == "Periodic") Grid->PeriodicBoundaryConditions(systemMaxwell, b_vector,'z', low_rank, high_rank);

                // Apply source
                for (size_t i = 0; i < config.sourceLocations.size(); ++i) {
                    const auto& loc = config.sourceLocations[i];
                    std::vector<Complex> value = config.sourceValues[i];
            
                    for(int j = 0; j<3; j++){
                        if (std::abs(value[j].real() - Constants::FREE.real()) > 1e-10 && std::abs(value[j].imag() - Constants::FREE.imag()) > 1e-10){
                            temp_idx = Grid->point_index(std::get<0>(loc), std::get<1>(loc), std::get<2>(loc)) + j;
                            if(temp_idx >= low_rank && temp_idx < high_rank) VecSetValue(b_vector, temp_idx , value[j], INSERT_VALUES);
                        }
                    }
                }

                // Apply excitations
                for (size_t i = 0; i < config.excitationLocations.size(); ++i) {
                    const auto& loc = config.excitationLocations[i];
                    temp_idx = Grid->point_index(std::get<0>(loc), std::get<1>(loc), std::get<2>(loc));
                    std::vector<Complex> value = config.excitationValues[i];
                            
                    for (int j = 0; j<3; j++){
                        if (std::abs(value[j].real() - Constants::FREE.real()) > 1e-10 && std::abs(value[j].imag() - Constants::FREE.imag()) > 1e-10){
                            Grid->DirichletBoundaryConditionsPoint(systemMaxwell, b_vector, value[j], temp_idx + j, low_rank, high_rank);
                        } 
                    }
                }

                Grid->AddMetalPointsInVolume(systemMaxwell, b_vector, config.metalLocations, config.y_reflector, low_rank, high_rank);

                MatAssemblyBegin(systemMaxwell, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(systemMaxwell, MAT_FINAL_ASSEMBLY);

                VecAssemblyBegin(b_vector);
                VecAssemblyEnd(b_vector);


                if(config.printMatrix) {
                    if(config.printProgress  && rank == 0) std::cout << "Printing Matrix" << std::endl;
                    printMatrixToFile(systemMaxwell, config.outputDirectory+"matrix.txt");
                    printVecToFile(b_vector, config.outputDirectory+"b_vector.txt");
                }

                MPI_Barrier(MPI_COMM_WORLD);

                if(print_times){
                    PetscTime(&endTime);
                    PetscPrintf(PETSC_COMM_WORLD, "Time taken on Applying Boundary Conditions: %.3f seconds\n", endTime - previous_time);
                    previous_time = endTime;
                }

                if(config.printProgress  && rank == 0) std::cout << "Solving Equations" << std::endl;

                MatrixSolver solver(num_variables, config.tolerance, config.maxIterations, config.scalar);
                if(config.solverMethod != "direct_LU" && config.initial_guess_file != "") solver.setInitialGuessFromFileBinary(harps_dir + config.initial_guess_file);
                
                solver.solve(systemMaxwell, b_vector, size, config.solverMethod, config.PC); // "sor", "amg"
                solver.printConvergenceInfo();

                Vec solution = solver.getSolution();

                if(print_times){
                    PetscTime(&endTime);
                    PetscPrintf(PETSC_COMM_WORLD, "Time taken on Solving System: %.3f seconds\n", endTime - previous_time);
                    previous_time = endTime;
                }

                if(config.printSolution){
                    printVecToFileBinary(solution, config.outputDirectory + "fields_solution.bin");
                    //printVecToFile(solution, config.outputDirectory + "fields_solution.txt");
                    
                    if(print_times){
                        PetscTime(&endTime);
                        PetscPrintf(PETSC_COMM_WORLD, "Time taken on Storing Solution: %.3f seconds\n", endTime - previous_time);
                        previous_time = endTime;
                    }
                }

                // Create a sequential vector to hold the complete solution on each process
                MPI_Barrier(MPI_COMM_WORLD);
                if(config.printProgress  && rank == 0) std::cout << "Creating Sequential Vector" << std::endl;

                Vec x_seq;
                VecCreate(PETSC_COMM_SELF, &x_seq);
                VecSetSizes(x_seq, PETSC_DECIDE, num_variables);
                VecSetFromOptions(x_seq);

                VecScatter scatter;
                VecScatterCreateToAll(solution, &scatter, &x_seq);

                VecScatterBegin(scatter, solution, x_seq, INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd(scatter, solution, x_seq, INSERT_VALUES, SCATTER_FORWARD);

                Complex* fields = new Complex[num_variables];

                // Convert scalar field to vector field for storing
                if(config.scalar != -1){
                    Complex* fields_scalar;
                    VecGetArray(x_seq, &fields_scalar);
                    for(int i = num_pts-1; i > -1; i--) fields[i*3+config.scalar] = fields_scalar[i];

                    // set the rest of the components to zero
                    for(int i = 0; i < num_pts; i++) {
                        fields[i*3 + (config.scalar+1)%3] = zero_C;
                        fields[i*3 + (config.scalar+2)%3] = zero_C;
                    }
                } else{
                    VecGetArray(x_seq, &fields);
                }

                if(config.printProgress  && rank == 0) std::cout << "Calculating Absorbed Power" << std::endl;
                
                double* absorbedPowerDensity;
                if(rank == 0){
                    absorbedPowerDensity = Grid->computeAbsorbedPowerDens(fields, real_conductivity, num_pts);  // Using Joule Heating
                    // absorbedPowerDensity = Grid->normalizeDens(absorbedPowerDensity, num_pts);                  // Normalize P_abs to 1
                    Grid->normalizeDens(absorbedPowerDensity, num_pts, total_abs_power); // without normalizing but still storing p_abs.txt
                    if(symmetric_harps) total_abs_power *= 2.0; // account for the other half of the domain

                    std::cout << "Total Joule Power: " << total_abs_power << std::endl;

                    // Converting from 2D into 1D radial and compacting into vectors to be sent to macros
                    double dr = 0.0004;
                    int num_bins = static_cast<int>(std::ceil(R_in / dr));

                    r_bins.resize(num_bins + 1);
                    for (int i = 0; i <= num_bins; ++i)  r_bins[i] = i * dr;

                    r_centers.resize(num_bins);
                    for (int i = 0; i < num_bins; ++i) r_centers[i] = (i + 0.5) * dr;  

                    std::vector<std::vector<double>> radial_power(N_macro, std::vector<double>(num_bins, 0.0));
                    std::vector<std::vector<double>> radial_ne(N_macro, std::vector<double>(num_bins, 0.0));
                    std::vector<std::vector<double>> radial_volume(N_macro, std::vector<double>(num_bins, 0.0));

                    for(int i = 0; i < num_pts; i++){    
                        double ne_temp;                    
                        x = Grid->getX( Grid->getXindex(i*3));
                        y = Grid->getY(Grid->getYindex(i*3));
                        //z = Grid->getZ(Grid->getZindex(i*3));
                        r = std::sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y));
                        if(r > R_in) continue;

                        theta = angle(x - center_x, y - center_y);

                        volume = Grid->getVolume(Grid->getXindex(i*3),Grid->getYindex(i*3),Grid->getZindex(i*3));
                        pabs = absorbedPowerDensity[i];
                        ne_temp = electron_density[i];

                        // Bin
                        int bin_idx = std::lower_bound(r_bins.begin(), r_bins.end(), r) - r_bins.begin() - 1;

                        if(theta >= 2*M_PI-angle_macro_2 || theta <= angle_macro_2){
                            radial_power[0][bin_idx] += pabs * volume;
                            radial_ne[0][bin_idx] += ne_temp * volume;
                            radial_volume[0][bin_idx] += volume;
                        }else{
                            for (int j = 1; j < N_macro; j++){
                                if(theta < (2*j+1)*angle_macro_2){
                                    radial_power[j][bin_idx] += pabs * volume;
                                    radial_ne[j][bin_idx] += ne_temp * volume;
                                    radial_volume[j][bin_idx] += volume;

                                    break;
                                }
                            }
                        }
                    }

                    std::vector<std::vector<double>> radial_power_density(N_macro, std::vector<double>(num_bins, 0.0));
                    std::vector<std::vector<double>> radial_ne_density(N_macro, std::vector<double>(num_bins, 0.0));
                    for (int i = 0; i < num_bins; ++i) {
                        for (int j = 0; j < N_macro; j++) {
                            if (radial_volume[j][i] > 0){
                                radial_power_density[j][i] = radial_power[j][i] / radial_volume[j][i];
                                radial_ne_density[j][i] = radial_ne[j][i] / radial_volume[j][i];
                            }                            
                        }
                    }

                    // Normalize p(r) 
                    double p_r_integral = 0.0;
                    for (int j = 0; j < N_macro; j ++){
                        for (int i = 0; i < num_bins; ++i) {
                            double r1 = r_bins[i];
                            double r2 = r_bins[i+1];
                            double r_mid = 0.5*(r1+r2);
                            double shell_area = r_mid*(r2-r1);  // annulus area
                            p_r_integral += radial_power_density[j][i] * shell_area;
                        }
                    }

                    // Pack into flat vector r_P_ne_out[i] = {r1[i], P1[i], r2[i], P2[i], ...}
                    for (int j = 0; j < N_macro; j++){
                        r_P_ne_out[j].clear();
                        r_P_ne_out[j].reserve(3 * num_bins);
                        for (int i = 0; i < num_bins; ++i) {
                            r_P_ne_out[j].push_back(r_centers[i]);
                            r_P_ne_out[j].push_back(N_macro*radial_power_density[j][i]/p_r_integral);
                            r_P_ne_out[j].push_back(radial_ne_density[j][i]);
                        }
                    }

                    if(config.storeResults) {
                        if(config.printProgress) std::cout << "Storing results" << std::endl;

                        // Helper function to check if a field should be exported
                        auto shouldExport = [](const std::string& fieldName, const std::vector<std::string>& outputFields) {
                            return std::find(outputFields.begin(), outputFields.end(), fieldName) != outputFields.end();
                        };            


                        if(shouldExport("PowerFlux",config.outputFields)){
                            double* P_flux = new double[3*num_pts];
                            double* P_residual = new double[num_pts];
                            double* P_flux_amplitude = new double[num_pts];
                            double* P_flux_x = new double[num_pts];
                            double* P_flux_y = new double[num_pts];
                            double* P_flux_z = new double[num_pts];

                            Grid->calculateFlux(fields, P_flux, angular_frequency, num_pts);
                            Grid->calculateResidual(P_flux, absorbedPowerDensity, P_residual, num_pts);
                            
                            Grid->calculateFieldAmplitudes(P_flux, P_flux_amplitude, num_pts);
                            Grid->getFieldComponents(P_flux, P_flux_x, P_flux_y, P_flux_z, num_pts);

                            exportField3D(P_flux_amplitude, config.outputDirectory + "power_flux_amp.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);    
                            exportField3D(P_residual, config.outputDirectory + "power_residual.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ , non_uniform_grid, Grid->x, Grid->y, Grid->z);

                            exportField3D(P_flux_x, config.outputDirectory + "power_flux_x.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z); 
                            exportField3D(P_flux_y, config.outputDirectory + "power_flux_y.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);  
                            exportField3D(P_flux_z, config.outputDirectory + "power_flux_z.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            
                            delete[] P_flux_amplitude; delete[] P_flux_x; delete[] P_flux_y; delete[] P_flux_z; delete[] P_flux; delete[] P_residual;
                        }
                
                        if(shouldExport("AbsorbedPower",config.outputFields)){
                            exportField3D(absorbedPowerDensity, config.outputDirectory + "absorbed_power_density.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);    
                        }
                        if(shouldExport("Inputs",config.outputFields)){
                            exportField3D(electron_density, config.outputDirectory + "electron_density.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ , non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(real_mobility, config.outputDirectory + "real_mobility.txt", config.n_x, config.n_y, config.n_z,
                                    config.lengthX, config.lengthY, config.lengthZ , non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(imag_mobility, config.outputDirectory + "imag_mobility.txt", config.n_x, config.n_y, config.n_z,
                                    config.lengthX, config.lengthY, config.lengthZ , non_uniform_grid, Grid->x, Grid->y, Grid->z);
                        }
                        if(shouldExport("FieldAmplitudes",config.outputFields)){
                            double* E_amplitude = new double[num_pts];
                            Grid->calculateFieldAmplitudes(fields, E_amplitude, num_pts);

                            exportField3D(E_amplitude, config.outputDirectory + "E_amplitude.txt", config.n_x, config.n_y, config.n_z, 
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            
                            delete[] E_amplitude;
                        }
                        if(shouldExport("FieldRealPart",config.outputFields)){
                            double* E_real = new double[num_pts];
                            double* E_imag = new double[num_pts];
                            Grid->calculateFieldRealPart(fields, E_real, num_pts);
                            Grid->calculateFieldImagPart(fields, E_imag, num_pts);

                            exportField3D(E_real, config.outputDirectory + "E_real.txt", config.n_x, config.n_y, config.n_z, 
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(E_imag, config.outputDirectory + "E_imag.txt", config.n_x, config.n_y, config.n_z, 
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                        
                            delete[] E_real;    delete[] E_imag;
                        }
                        if(shouldExport("FieldComponents",config.outputFields)){
                            double* Ex_amp = new double[num_pts];
                            double* Ey_amp = new double[num_pts];
                            double* Ez_amp = new double[num_pts];
                            Grid->getFieldComponents(fields, Ex_amp, Ey_amp, Ez_amp, num_pts);

                            exportField3D(Ex_amp, config.outputDirectory + "Ex_amp.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(Ey_amp, config.outputDirectory + "Ey_amp.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(Ez_amp, config.outputDirectory + "Ez_amp.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);

                            delete[] Ex_amp; delete[] Ey_amp; delete[] Ez_amp;
                        }
                        if(shouldExport("F_Grad_Components",config.outputFields)){
                            double* f_grad_cond_x = new double[num_pts];
                            double* f_grad_cond_y = new double[num_pts];
                            double* f_grad_cond_z = new double[num_pts];
                            Grid->getFieldComponents(f_grad_cond,f_grad_cond_x, f_grad_cond_y, f_grad_cond_z, num_pts);

                            exportField3D(f_grad_cond_x, config.outputDirectory + "f_grad_cond_x.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(f_grad_cond_y, config.outputDirectory + "f_grad_cond_y.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                            exportField3D(f_grad_cond_z, config.outputDirectory + "f_grad_cond_z.txt", config.n_x, config.n_y, config.n_z,
                                config.lengthX, config.lengthY, config.lengthZ, non_uniform_grid, Grid->x, Grid->y, Grid->z);
                        
                            delete[] f_grad_cond_x; delete[] f_grad_cond_y; delete[] f_grad_cond_z;
                        }
                    }

                    if(print_times){
                        PetscTime(&endTime);
                        PetscPrintf(PETSC_COMM_WORLD, "Time taken on Storing Results and Calculating Pabs: %.3f seconds\n", endTime - previous_time);
                        previous_time = endTime;
                    }
                    if(config.printProgress) std::cout << "End of Program" << std::endl;

                    delete absorbedPowerDensity;
                }
                PetscTime(&endTime);
                PetscPrintf(PETSC_COMM_WORLD, "Time taken: %.3f seconds\n", endTime - startTime);

                delete[] electron_density; delete[] real_mobility; delete[] imag_mobility; delete[] fields;  delete[] f_grad_cond;
                delete[] real_conductivity; delete[] complex_conductivity; delete[] complex_permittivity; delete[] real_permittivity;

                VecScatterDestroy(&scatter);
                VecDestroy(&x_seq);
                VecDestroy(&solution);
                VecDestroy(&b_vector);
                MatDestroy(&systemMaxwell);
                
            } catch (const std::exception& error_config) {
                std::cerr << "Error: " << error_config.what() << std::endl;
            }
        }

        // O_F
        if(rank == 0){
            for(int i = 0; i < N_macro; ++i) {
                if(do_harps){
                    auto harps_result = Data::grid(r_P_ne_out[i].data(), {r_P_ne_out[i].size()}, {"rPne"});
                    instance->send( "end_harps_" + std::to_string(i), Message(total_abs_power, harps_result));
                }else{
                    instance->send("end_harps_" + std::to_string(i), Message(total_abs_power, Data(flag_sync)));
                }
            }
        }
    }

    PetscFinalize();
}

int main(int argc, char* argv[]) {
    harps_main(argc, argv);
    return EXIT_SUCCESS;
}