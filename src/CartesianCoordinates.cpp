#include "CartesianCoordinates.h"
#include <fstream>
#include <algorithm>

//Variable in the system of equations: (E_x, E_y, E_z) for each point (i,j,k)
//Point index: 3*(i*n_y*n_z + j*n_z + k)

Mat CartesianCoordinateSystem::createMaxwellEquationMatrix(Complex* f_grad_cond, Complex* complex_permittivity, double vacuum_wave_number, double waveguide_number, double yCenter, int size){
    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n_variables, n_variables);  
    MatSetFromOptions(A);
    
    if(size > 1) MatMPIAIJSetPreallocation(A, 15, NULL, 15, NULL); 
    else MatSeqAIJSetPreallocation(A, 15, NULL);
    
    MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
    MatSetBlockSize(A, 3);
    MatSetUp(A);

    // Get the rows owned by this rank so matrix is not calculated multiple times
    PetscInt rstart, rend;
    MatGetOwnershipRange(A, &rstart, &rend);

    PetscInt row;
    Complex h_inv_sq, main_diag_x, main_diag_y, main_diag_z;

    // Applying PML coordinate stretching
    dx_s.clear(); dy_s.clear(); dz_s.clear();
    for (int i = 0; i < n_x-1; i++) dx_s.push_back(dx[i] * s_profile_x[i]);
    for (int j = 0; j < n_y-1; j++) dy_s.push_back(dy[j] * s_profile_y[j]);
    for (int k = 0; k < n_z-1; k++) dz_s.push_back(dz[k] * s_profile_z[k]);

    // Reusable vectors for batched row insertion
    std::vector<PetscInt> cols;
    std::vector<Complex> vals;
    cols.reserve(15);
    vals.reserve(15);


    // Fill A with system from Maxwell ∇²E + ∇(f·E) + gE = 0

    // E is a vector, ∇(f·E) is gradient of inner product
    // f = ∇(cond)/(i ω ε_0 ε_r_complex)
    // g = k_0 ^2 * ε_r_complex
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            for (int k = 0; k < n_z; k++) {
                
                row = point_index(i, j, k);

                // CRITICAL GUARD: Only the process owning these rows does the work
                if (row >= rstart && row < rend) {
                    h_inv_sq = 0; main_diag_x = 0; main_diag_y = 0; main_diag_z = 0;

                    if(i > 0 && i < n_x-1){
                        h_inv_sq += 1./(dx_s[i]*dx_s[i-1]);     // ∇²E_xx
                        main_diag_x = (1./dx_s[i-1] - 1./dx_s[i])*f_grad_cond[row]; // ∇(f·E)_xx
                    }
                    if(j > 0 && j < n_y-1){
                        h_inv_sq += 1./(dy_s[j]*dy_s[j-1]);     // ∇²E_yy
                        main_diag_y = (1./dy_s[j-1] - 1./dy_s[j])*f_grad_cond[row+1]; // ∇(f·E)_yy
                    }
                    if(k > 0 && k < n_z-1){
                        h_inv_sq += 1./(dz_s[k]*dz_s[k-1]);     // ∇²E_zy
                        main_diag_z = (1./dz_s[k-1] - 1./dz_s[k])*f_grad_cond[row+2]; // ∇(f·E)_zz
                    }
                    
                    main_diag_x += -2*h_inv_sq + vacuum_wave_number*vacuum_wave_number*complex_permittivity[row/3]; 
                    main_diag_y += -2*h_inv_sq + vacuum_wave_number*vacuum_wave_number*complex_permittivity[row/3];
                    main_diag_z += -2*h_inv_sq + vacuum_wave_number*vacuum_wave_number*complex_permittivity[row/3];

                    if(n_x == 1){
                        if(abs(y[j] - yCenter) > 0.013) main_diag_z -= waveguide_number*waveguide_number;  // When n_x = 1, we still consider TE_10 influence in the original wave inside the waveguide
                    }

                    // Off main diagonal terms ∇²E + ∇(f·E)
                    // On E_x
                    cols.clear(); vals.clear();
                    if (i > 0 && i < n_x - 1) {
                        // ∇²E_xx + ∇(f·E)_xx
                        cols.push_back(point_index(i - 1, j, k)); 
                        vals.push_back(2. / (dx_s[i - 1] * (dx_s[i - 1] + dx_s[i])) + (1. / (dx_s[i - 1] + dx_s[i]) - 1. / dx_s[i - 1]) * f_grad_cond[point_index(i - 1, j, k)]);
                        
                        cols.push_back(point_index(i + 1, j, k)); 
                        vals.push_back(2. / (dx_s[i] * (dx_s[i - 1] + dx_s[i])) + (1. / dx_s[i] - 1. / (dx_s[i - 1] + dx_s[i])) * f_grad_cond[point_index(i + 1, j, k)]);

                        // ∇(f·E)_yx
                        cols.push_back(point_index(i, j, k) + 1);     vals.push_back((1. / dx_s[i - 1] - 1. / dx_s[i]) * f_grad_cond[row + 1]);
                        cols.push_back(point_index(i - 1, j, k) + 1); vals.push_back((1. / (dx_s[i - 1] + dx_s[i]) - 1. / dx_s[i - 1]) * f_grad_cond[point_index(i - 1, j, k) + 1]);
                        cols.push_back(point_index(i + 1, j, k) + 1); vals.push_back((1. / dx_s[i] - 1. / (dx_s[i - 1] + dx_s[i])) * f_grad_cond[point_index(i + 1, j, k) + 1]);

                        // ∇(f·E)_zx
                        cols.push_back(point_index(i, j, k) + 2);     vals.push_back((1. / dx_s[i - 1] - 1. / dx_s[i]) * f_grad_cond[row + 2]);
                        cols.push_back(point_index(i - 1, j, k) + 2); vals.push_back((1. / (dx_s[i - 1] + dx_s[i]) - 1. / dx_s[i - 1]) * f_grad_cond[point_index(i - 1, j, k) + 2]);
                        cols.push_back(point_index(i + 1, j, k) + 2); vals.push_back((1. / dx_s[i] - 1. / (dx_s[i - 1] + dx_s[i])) * f_grad_cond[point_index(i + 1, j, k) + 2]);
                    }
                    if (j > 0 && j < n_y - 1) {
                        // ∇²E_xy
                        cols.push_back(point_index(i, j - 1, k)); vals.push_back(2. / (dy_s[j - 1] * (dy_s[j - 1] + dy_s[j]))); 
                        cols.push_back(point_index(i, j + 1, k)); vals.push_back(2. / (dy_s[j] * (dy_s[j - 1] + dy_s[j])));
                    }
                    if (k > 0 && k < n_z - 1) {
                        // ∇²E_xz
                        cols.push_back(point_index(i, j, k - 1)); vals.push_back(2. / (dz_s[k - 1] * (dz_s[k - 1] + dz_s[k]))); 
                        cols.push_back(point_index(i, j, k + 1)); vals.push_back(2. / (dz_s[k] * (dz_s[k - 1] + dz_s[k])));
                    }
                    // Single batched call for E_x row
                    cols.push_back(row); vals.push_back(main_diag_x);

                    MatSetValues(A, 1, &row, cols.size(), cols.data(), vals.data(), INSERT_VALUES);

                    // On E_y
                    cols.clear(); vals.clear();
                    PetscInt row_y = row + 1;
                    if (i > 0 && i < n_x - 1) {
                        // ∇²E_yx
                        cols.push_back(point_index(i - 1, j, k) + 1); vals.push_back(2. / (dx_s[i - 1] * (dx_s[i - 1] + dx_s[i]))); 
                        cols.push_back(point_index(i + 1, j, k) + 1); vals.push_back(2. / (dx_s[i] * (dx_s[i - 1] + dx_s[i])));
                    }
                    if (j > 0 && j < n_y - 1) {
                        // ∇²E_yy + ∇(f·E)_yy
                        cols.push_back(point_index(i, j - 1, k) + 1); 
                        vals.push_back(2. / (dy_s[j - 1] * (dy_s[j - 1] + dy_s[j])) + (1. / (dy_s[j - 1] + dy_s[j]) - 1. / dy_s[j - 1]) * f_grad_cond[point_index(i, j - 1, k) + 1]);
                        
                        cols.push_back(point_index(i, j + 1, k) + 1); 
                        vals.push_back(2. / (dy_s[j] * (dy_s[j - 1] + dy_s[j])) + (1. / dy_s[j] - 1. / (dy_s[j - 1] + dy_s[j])) * f_grad_cond[point_index(i, j + 1, k) + 1]);

                        // ∇(f·E)_xy
                        cols.push_back(point_index(i, j, k) + 0);     vals.push_back((1. / dy_s[j - 1] - 1. / dy_s[j]) * f_grad_cond[row]);
                        cols.push_back(point_index(i, j - 1, k) + 0); vals.push_back((1. / (dy_s[j - 1] + dy_s[j]) - 1. / dy_s[j - 1]) * f_grad_cond[point_index(i, j - 1, k)]);
                        cols.push_back(point_index(i, j + 1, k) + 0); vals.push_back((1. / dy_s[j] - 1. / (dy_s[j - 1] + dy_s[j])) * f_grad_cond[point_index(i, j + 1, k)]);

                        // ∇(f·E)_zy
                        cols.push_back(point_index(i, j, k) + 2);     vals.push_back((1. / dy_s[j - 1] - 1. / dy_s[j]) * f_grad_cond[row + 2]);
                        cols.push_back(point_index(i, j - 1, k) + 2); vals.push_back((1. / (dy_s[j - 1] + dy_s[j]) - 1. / dy_s[j - 1]) * f_grad_cond[point_index(i, j - 1, k) + 2]);
                        cols.push_back(point_index(i, j + 1, k) + 2); vals.push_back((1. / dy_s[j] - 1. / (dy_s[j - 1] + dy_s[j])) * f_grad_cond[point_index(i, j + 1, k) + 2]);
                    }
                    if (k > 0 && k < n_z - 1) {
                        // ∇²E_yz
                        cols.push_back(point_index(i, j, k - 1) + 1); vals.push_back(2. / (dz_s[k - 1] * (dz_s[k - 1] + dz_s[k]))); 
                        cols.push_back(point_index(i, j, k + 1) + 1); vals.push_back(2. / (dz_s[k] * (dz_s[k - 1] + dz_s[k])));
                    }
                    cols.push_back(row_y); vals.push_back(main_diag_y);

                    MatSetValues(A, 1, &row_y, cols.size(), cols.data(), vals.data(), INSERT_VALUES);

                    // On E_z
                    cols.clear(); vals.clear();
                    PetscInt row_z = row + 2;
                    if (i > 0 && i < n_x - 1) {
                        // ∇²E_zx
                        cols.push_back(point_index(i - 1, j, k) + 2); vals.push_back(2. / (dx_s[i - 1] * (dx_s[i - 1] + dx_s[i]))); 
                        cols.push_back(point_index(i + 1, j, k) + 2); vals.push_back(2. / (dx_s[i] * (dx_s[i - 1] + dx_s[i])));
                    }
                    if (j > 0 && j < n_y - 1) {
                        // ∇²E_zy
                        cols.push_back(point_index(i, j - 1, k) + 2); vals.push_back(2. / (dy_s[j - 1] * (dy_s[j - 1] + dy_s[j]))); 
                        cols.push_back(point_index(i, j + 1, k) + 2); vals.push_back(2. / (dy_s[j] * (dy_s[j - 1] + dy_s[j])));
                    } 
                    if (k > 0 && k < n_z - 1) {
                        // ∇²E_zz + ∇(f·E)_zz
                        cols.push_back(point_index(i, j, k - 1) + 2); 
                        vals.push_back(2. / (dz_s[k - 1] * (dz_s[k - 1] + dz_s[k])) + (1. / (dz_s[k - 1] + dz_s[k]) - 1. / dz_s[k - 1]) * f_grad_cond[point_index(i, j, k - 1) + 2]);
                        
                        cols.push_back(point_index(i, j, k + 1) + 2); 
                        vals.push_back(2. / (dz_s[k] * (dz_s[k - 1] + dz_s[k])) + (1. / dz_s[k] - 1. / (dz_s[k - 1] + dz_s[k])) * f_grad_cond[point_index(i, j, k + 1) + 2]);

                        // ∇(f·E)_xz
                        cols.push_back(point_index(i, j, k) + 0);     vals.push_back((1. / dz_s[k - 1] - 1. / dz_s[k]) * f_grad_cond[row + 0]);
                        cols.push_back(point_index(i, j, k - 1) + 0); vals.push_back((1. / (dz_s[k - 1] + dz_s[k]) - 1. / dz_s[k - 1]) * f_grad_cond[point_index(i, j, k - 1)]);
                        cols.push_back(point_index(i, j, k + 1) + 0); vals.push_back((1. / dz_s[k] - 1. / (dz_s[k - 1] + dz_s[k])) * f_grad_cond[point_index(i, j, k + 1)]);

                        // ∇(f·E)_yz
                        cols.push_back(point_index(i, j, k) + 1);     vals.push_back((1. / dz_s[k - 1] - 1. / dz_s[k]) * f_grad_cond[row + 1]);
                        cols.push_back(point_index(i, j, k - 1) + 1); vals.push_back((1. / (dz_s[k - 1] + dz_s[k]) - 1. / dz_s[k - 1]) * f_grad_cond[point_index(i, j, k - 1) + 1]);
                        cols.push_back(point_index(i, j, k + 1) + 1); vals.push_back((1. / dz_s[k] - 1. / (dz_s[k - 1] + dz_s[k])) * f_grad_cond[point_index(i, j, k + 1) + 1]);
                    }
                    cols.push_back(row_z); vals.push_back(main_diag_z);
                
                    MatSetValues(A, 1, &row_z, cols.size(), cols.data(), vals.data(), INSERT_VALUES);
                }
            }
        }
    }
    
    return A;
}


std::vector<Complex> CartesianCoordinateSystem::calculateCondGradFunction(Complex* complex_conductivity, Complex* complex_permittivity, double angular_frequency){
    std::vector<Complex> f_grad_cond(n_points*3, Complex(0., 0.));

    Complex grad_x, grad_y, grad_z;
    int component_idx;

    for (int i = 0; i < n_x; ++i) {
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                component_idx = point_index(i,j,k);
                grad_x = 0; grad_y = 0; grad_z = 0;

                // Calculating ∇(cond) with central differences in non-uniform grid
                if(i > 0 && i < n_x-1){
                    grad_x = (1./dx[i] - 1./(dx[i-1]+dx[i])) * complex_conductivity[point_index(i+1,j,k)/3] 
                           + (1./dx[i-1] - 1./dx[i]) * complex_conductivity[point_index(i,j,k)/3]
                           + (1./(dx[i-1]+dx[i]) - 1./dx[i-1]) * complex_conductivity[point_index(i-1,j,k)/3];
                }
                if(j > 0 && j < n_y-1){
                    grad_y = (1./dy[j] - 1./(dy[j-1]+dy[j])) * complex_conductivity[point_index(i,j+1,k)/3] 
                           + (1./dy[j-1] - 1./dy[j]) * complex_conductivity[point_index(i,j,k)/3]
                           + (1./(dy[j-1]+dy[j]) - 1./dy[j-1]) * complex_conductivity[point_index(i,j-1,k)/3];
                }
                if(k > 0 && k < n_z-1){
                    grad_z = (1./dz[k] - 1./(dz[k-1]+dz[k])) * complex_conductivity[point_index(i,j,k+1)/3] 
                           + (1./dz[k-1] - 1./dz[k]) * complex_conductivity[point_index(i,j,k)/3]
                           + (1./(dz[k-1]+dz[k]) - 1./dz[k-1]) * complex_conductivity[point_index(i,j,k-1)/3];
                }

                // f = ∇(cond) / (i ω ε_0 ε_r)
                Complex denominator = Complex(0,1) * angular_frequency * Constants::EPSILON_0 * complex_permittivity[point_index(i,j,k)/3];

                f_grad_cond[component_idx]     = grad_x / denominator;  // fx
                f_grad_cond[component_idx + 1] = grad_y / denominator;  // fy
                f_grad_cond[component_idx + 2] = grad_z / denominator;  // fz
            }
        }
    }

    return f_grad_cond;
}


void CartesianCoordinateSystem::DirichletBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction, bool upper_boundary, Complex value, PetscInt low_rank, PetscInt high_rank){
    int eq_point_index;

    int l_array[3];
    int l_size;
    for (int i = 0; i < 3; ++i) l_array[i] = 0;

    if (field_str == "all"){
        l_size = 3;
        for (int i = 0; i < l_size; ++i) l_array[i] = i;
    }else if(field_str == "E_x"){
        l_size = 1;
        l_array[0] = 0;
    }else if(field_str == "E_y"){
        l_size = 1;
        l_array[0] = 1;
    }else if(field_str == "E_z"){
        l_size = 1;
        l_array[0] = 2;
    }else{
        std::cout << "Invalid field string for Dirichlet Boundary Conditions (all, E_x, E_y, E_z)" << std::endl;
    }

    if(direction == 'x'){
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                if(upper_boundary) eq_point_index = point_index(n_x-1,j,k);
                else               eq_point_index = point_index(0,j,k);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    DirichletBoundaryConditionsPoint(matrix, b_array, value, eq_point_index+l_array[l_index], low_rank, high_rank);
                }
            }
        }
    }else if(direction == 'y'){
        for (int i = 0; i < n_x; ++i) {
            for (int k = 0; k < n_z; ++k) {
                if(upper_boundary) eq_point_index = point_index(i,n_y-1,k);
                else               eq_point_index = point_index(i,0,k);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    DirichletBoundaryConditionsPoint(matrix, b_array, value, eq_point_index+l_array[l_index], low_rank, high_rank);
                }
            }
        }
    }else if(direction == 'z'){
        for (int i = 0; i < n_x; ++i) {
            for (int j = 0; j < n_y; ++j) {
                if(upper_boundary) eq_point_index = point_index(i,j,n_z-1);
                else               eq_point_index = point_index(i,j,0);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    DirichletBoundaryConditionsPoint(matrix, b_array, value, eq_point_index+l_array[l_index], low_rank, high_rank);
                }
            }
        }
    }else{
        std::cout << "Invalid direction for Dirichlet Boundary Conditions (x,y,z)" << std::endl;
    }
}


void CartesianCoordinateSystem::DirichletBoundaryConditionsPoint(Mat& matrix, Vec& b_array, Complex value, int idx, PetscInt low_rank, PetscInt high_rank) {
    PetscInt row = idx;

    if (row >= low_rank && row < high_rank) {
        SetRowZeroManually(matrix, row, n_variables, 3);

        MatSetValue(matrix, row, row, 1.0, INSERT_VALUES);
        VecSetValue(b_array, row, value, INSERT_VALUES);
    }
}


void CartesianCoordinateSystem::NeumannBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction, bool upper_boundary, Complex value, PetscInt low_rank, PetscInt high_rank){
    int eq_point_index;

    int l_array[3];
    int l_size;
    for (int i = 0; i < 3; ++i) l_array[i] = 0;

    if (field_str == "all"){
        l_size = 3;
        for (int i = 0; i < l_size; ++i) l_array[i] = i;
    }else if(field_str == "E_x"){
        l_size = 1;
        l_array[0] = 0;
    }else if(field_str == "E_y"){
        l_size = 1;
        l_array[0] = 1;
    }else if(field_str == "E_z"){
        l_size = 1;
        l_array[0] = 2;
    }else{
        std::cout << "Invalid field string for Neumann Boundary Conditions (all, E_x, E_y, E_z)" << std::endl;
    }

    if(direction == 'x'){
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                if(upper_boundary) eq_point_index = point_index(n_x-1,j,k);
                else               eq_point_index = point_index(0,j,k);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    PetscInt row = eq_point_index + l_array[l_index];
                    if(row >= low_rank && row < high_rank) SetRowZeroManually(matrix, row, n_variables, 3);

                    if(upper_boundary && row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, 1., INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(n_x-2,j,k)+l_array[l_index], -1., INSERT_VALUES);
                    }else if(row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, -1., INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(1,j,k)+l_array[l_index], 1., INSERT_VALUES);
                    }
                    if(upper_boundary && row >= low_rank && row < high_rank )   VecSetValue(b_array, row, value*dx[n_x-2], INSERT_VALUES);
                    else if(row >= low_rank && row < high_rank)                 VecSetValue(b_array, row, value*dx[0], INSERT_VALUES);
                }
            }
        }
    }else if(direction == 'y'){
        for (int i = 0; i < n_x; ++i) {
            for (int k = 0; k < n_z; ++k) {
                if(upper_boundary) eq_point_index = point_index(i,n_y-1,k);
                else               eq_point_index = point_index(i,0,k);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    PetscInt row = eq_point_index + l_array[l_index];
                    if(row >= low_rank && row < high_rank) SetRowZeroManually(matrix, row, n_variables, 3);

                    if(upper_boundary && row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, 1., INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,n_y-2,k)+l_array[l_index], -1., INSERT_VALUES);
                    }else if(row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, -1., INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,1,k)+l_array[l_index], 1., INSERT_VALUES);
                    }
                    if(upper_boundary && row >= low_rank && row < high_rank)  VecSetValue(b_array, row, value*dy[n_y-2], INSERT_VALUES);
                    else if(row >= low_rank && row < high_rank)               VecSetValue(b_array, row, value*dy[0], INSERT_VALUES);
                }
            }
        }
    }else if(direction == 'z'){
        for (int i = 0; i < n_x; ++i) {
            for (int j = 0; j < n_y; ++j) {
                if(upper_boundary) eq_point_index = point_index(i,j,n_z-1);
                else               eq_point_index = point_index(i,j,0);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    PetscInt row = eq_point_index + l_array[l_index];
                    if(row >= low_rank && row < high_rank) SetRowZeroManually(matrix, row, n_variables, 3);

                    if(upper_boundary && row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, 1., INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,j,n_z-2)+l_array[l_index], -1., INSERT_VALUES);
                    }else if (row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, -1., INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,j,1)+l_array[l_index], 1., INSERT_VALUES);
                    }
                    if(upper_boundary && row >= low_rank && row < high_rank)    VecSetValue(b_array, row, value*dz[n_z-2], INSERT_VALUES);
                    else if(row >= low_rank && row < high_rank)                 VecSetValue(b_array, row, value*dz[0], INSERT_VALUES);
                }
            }
        }
    }else{
        std::cout << "Invalid direction for Neumann Boundary Conditions (x,y,z)" << std::endl;
    }
}


void CartesianCoordinateSystem::LocalNeumannCondition(Mat& matrix, Vec& b_array, Complex value, int idx_a, int idx_b, char direction, PetscInt low_rank, PetscInt high_rank){
    PetscInt row = idx_a;
    int delta;

    if(direction == 'x')        delta = dx[getXindex(idx_a)];
    else if(direction == 'y')   delta = dy[getYindex(idx_a)];
    else if(direction == 'z')   delta = dz[getZindex(idx_a)];
    else throw std::runtime_error("Invalid direction for Neumann Local Boundary Condition (x,y,z)");

    if(row >= low_rank && row < high_rank){
        SetRowZeroManually(matrix, row, n_variables, 3);
        MatSetValue(matrix, row, row, 1.0, INSERT_VALUES);
        MatSetValue(matrix, row, idx_b, -1.0, INSERT_VALUES);

        VecSetValue(b_array, row, value * delta, INSERT_VALUES);
    }
}


void CartesianCoordinateSystem::RobinBoundaryConditions(Mat& matrix, Vec& b_array, std::string field_str, char direction, bool upper_boundary,
                                std::vector<Complex> wave_value, double vacuum_wave_number, double cutoff_wave_number, PetscInt low_rank, PetscInt high_rank){
    int eq_point_index;

    int l_array[3];
    int l_size;

    for (int i = 0; i < 3; ++i) l_array[i] = 0;

    if (field_str == "all"){
        l_size = 3;
        for (int i = 0; i < l_size; ++i) l_array[i] = i;
    }else if(field_str == "E_x"){
        l_size = 1;
        l_array[0] = 0;
    }else if(field_str == "E_y"){
        l_size = 1;
        l_array[0] = 1;
    }else if(field_str == "E_z"){
        l_size = 1;
        l_array[0] = 2;
    }else{
        std::cout << "Invalid field string for Robin Boundary Conditions (all, E_x, E_y, E_z)" << std::endl;
    }

    double wave_number = sqrt(vacuum_wave_number*vacuum_wave_number - cutoff_wave_number*cutoff_wave_number); // beta for TE mode on waveguide
    //std::complex<double> wave_number = std::complex<double>(29.257909502765205,-2.538889203869415); // beta complex for TE mode on waveguide

    if(direction == 'x'){
        if((int) wave_value.size() != n_y*n_z) {
            std::cerr << "Error: Robin BC values size does not match grid dimensions. " << wave_value.size() << std::endl;
            return;
        }
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                if(upper_boundary) eq_point_index = point_index(n_x-1,j,k);
                else               eq_point_index = point_index(0,j,k);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    PetscInt row = eq_point_index + l_array[l_index];

                    if(((j == 0 || j == n_y-1 ) && n_x > 1) || ((k == 0 || k == n_z-1 ) && n_z > 1)){
                        //DirichletBoundaryConditionsPoint(matrix, b_array, wave_value[n_z*j+k], row, low_rank, high_rank);
                        continue;
                    }

                    if(row >= low_rank && row < high_rank) SetRowZeroManually(matrix, row, n_variables, 3);

                    if(upper_boundary && row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, 1./dx[n_x-2], INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(n_x-2,j,k)+l_array[l_index], Complex(0,1)*wave_number-1./dx[n_x-2], INSERT_VALUES);
                    }else if (row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, -Complex(0,1)*wave_number-1./dx[0], INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(1,j,k)+l_array[l_index], 1./dx[0], INSERT_VALUES);
                    }
                    if(upper_boundary && row >= low_rank && row < high_rank )   VecSetValue(b_array, row, 2*Complex(0,1)*wave_number*wave_value[n_z*j+k], INSERT_VALUES);
                    else if(row >= low_rank && row < high_rank)                 VecSetValue(b_array, row, -2*Complex(0,1)*wave_number*wave_value[n_z*j+k], INSERT_VALUES);
                }
            }
        }
    }else if(direction == 'y'){
        if((int) wave_value.size() != n_x*n_z) {
            std::cerr << "Error: Robin BC values size does not match grid dimensions. " << wave_value.size() << std::endl;
            return;
        }
        for (int i = 0; i < n_x; ++i) {
            for (int k = 0; k < n_z; ++k) {
                if(upper_boundary) eq_point_index = point_index(i,n_y-1,k);
                else               eq_point_index = point_index(i,0,k);

                for(int l_index = 0; l_index < l_size; l_index++){
                    PetscInt row = eq_point_index + l_array[l_index];

                    if(((i == 0 || i == n_x-1 ) && n_x > 1) || ((k == 0 || k == n_z-1 ) && n_z > 1)){
                        //DirichletBoundaryConditionsPoint(matrix, b_array, wave_value[n_z*i+k], row, low_rank, high_rank);
                        continue;
                    }

                    if(row >= low_rank && row < high_rank) SetRowZeroManually(matrix, row, n_variables, 3);

                    if(upper_boundary && row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, 1./dy[n_y-2], INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,n_y-2,k)+l_array[l_index], Complex(0,1)*wave_number-1./dy[n_y-2], INSERT_VALUES);
                    }else if(row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, -Complex(0,1)*wave_number-1./dy[0], INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,1,k)+l_array[l_index], 1./dy[0], INSERT_VALUES);
                    }
                    if(upper_boundary && row >= low_rank && row < high_rank) VecSetValue(b_array, row, 2*Complex(0,1)*wave_number*wave_value[n_z*i+k], INSERT_VALUES);
                    else if(row >= low_rank && row < high_rank)              VecSetValue(b_array, row, -2*Complex(0,1)*wave_number*wave_value[n_z*i+k], INSERT_VALUES);
                }
            }
        }
    }else if(direction == 'z'){
        if((int) wave_value.size() != n_x*n_y) {
            std::cerr << "Error: Robin BC values size does not match grid dimensions. " << wave_value.size() << std::endl;
            return;
        }
        for (int i = 0; i < n_x; ++i) {
            for (int j = 0; j < n_y; ++j) {
                if(upper_boundary) eq_point_index = point_index(i,j,n_z-1);
                else               eq_point_index = point_index(i,j,0);
                
                for(int l_index = 0; l_index < l_size; l_index++){
                    PetscInt row = eq_point_index + l_array[l_index];

                    if(((i == 0 || i == n_x-1 ) && n_x > 1) || ((j == 0 || j == n_y-1 ) && n_y > 1)){
                        //DirichletBoundaryConditionsPoint(matrix, b_array, wave_value[n_y*i+j], row, low_rank, high_rank);
                        continue;
                    }
                    if(row >= low_rank && row < high_rank) SetRowZeroManually(matrix, row, n_variables, 3);

                    if(upper_boundary && row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, 1./dz[n_z-2], INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,j,n_z-2)+l_array[l_index], Complex(0,1)*wave_number-1./dz[n_z-2], INSERT_VALUES);
                    }else if(row >= low_rank && row < high_rank){
                        MatSetValue(matrix, row, row, -Complex(0,1)*wave_number-1./dz[0], INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i,j,1)+l_array[l_index], 1./dz[0], INSERT_VALUES);
                    }
                    if(upper_boundary && row >= low_rank && row < high_rank)    VecSetValue(b_array, row, 2*Complex(0,1)*wave_number*wave_value[n_y*i+j], INSERT_VALUES);
                    else if(row >= low_rank && row < high_rank)                 VecSetValue(b_array, row, -2*Complex(0,1)*wave_number*wave_value[n_y*i+j], INSERT_VALUES);
                }
            }
        }
    }else{
        std::cout << "Invalid direction for Robin Boundary Conditions (x,y,z)" << std::endl;
    }
}


void CartesianCoordinateSystem::PeriodicBoundaryConditions(Mat& matrix, Vec& b_array, char direction, PetscInt low_rank, PetscInt high_rank){
    int eq_point_index;
    PetscInt row;
    
    if(direction == 'x'){        
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                eq_point_index = point_index(n_x-1, j, k);
                for(int l_index = 0; l_index < 3; l_index++){
                    row = eq_point_index + l_index;
                    if(row >= low_rank && row < high_rank) {
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, 1, INSERT_VALUES);
                    }
                }
            }
        }
        
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                eq_point_index = point_index(n_x-1, j, k);
                
                for(int l_index = 0; l_index < 3; l_index++) {
                    row = eq_point_index + l_index;
                    if(row >= low_rank && row < high_rank) {
                        MatSetValue(matrix, row, point_index(0, j, k) + l_index, -1., INSERT_VALUES);
                        VecSetValue(b_array, row, 0., INSERT_VALUES);
                    }
                }
            }
        }
    }else if(direction == 'y'){        
        for (int i = 0; i < n_x; ++i) {
            for (int k = 0; k < n_z; ++k) {
                eq_point_index = point_index(i,n_y-1, k);
                for(int l_index = 0; l_index < 3; l_index++){
                    row = eq_point_index + l_index;
                    if(row >= low_rank && row < high_rank) {
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, 1, INSERT_VALUES);
                    }
                }
            }
        }

        for (int i = 0; i < n_x; ++i) {
            for (int k = 0; k < n_z; ++k) {
                eq_point_index = point_index(i, n_y-1, k);
                
                for(int l_index = 0; l_index < 3; l_index++) {
                    row = eq_point_index + l_index;
                    if(row >= low_rank && row < high_rank) {
                        MatSetValue(matrix, row, point_index(i, 0, k) + l_index, -1., INSERT_VALUES);
                        VecSetValue(b_array, row, 0., INSERT_VALUES);
                    }
                }
            }
        }
    }else if(direction == 'z'){
        for (int i = 0; i < n_x; ++i) {
            for (int j = 0; j < n_y; ++j) {
                eq_point_index = point_index(i,j,n_z-1);
                for(int l_index = 0; l_index < 3; l_index++){
                    row = eq_point_index + l_index;
                    if(row >= low_rank && row < high_rank) {
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, 1, INSERT_VALUES);
                    }
                } 
            }
        }
        
        for (int i = 0; i < n_x; ++i) {
            for (int j = 0; j < n_y; ++j) {
                eq_point_index = point_index(i, j, n_z-1);
                
                for(int l_index = 0; l_index < 3; l_index++) {
                    row = eq_point_index + l_index;
                    if(row >= low_rank && row < high_rank) {
                        MatSetValue(matrix, row, point_index(i, j, 0) + l_index, -1., INSERT_VALUES);
                        VecSetValue(b_array, row, 0., INSERT_VALUES);
                    }
                }
            }
        }
    }else{
        std::cout << "Invalid direction for Periodic Boundary Conditions. Please try another one (x, y or z)" << std::endl;
    }
}


void CartesianCoordinateSystem::createPMLProfile(char direction, bool upper_boundary, int n_pml, double omega, int m, double sigma_0) {
    double x, sigma, sigma_factor, wavelenght;
    Complex s;

    wavelenght = 2*M_PI*Constants::C_LIGHT/omega;

    //sigma_factor = (m + 1) * 0.5 * 16.118 * Constants::C_LIGHT * wavelenght; // Perfectly Matched Layer (PML) for Computational Electromagnetics (2007) - Jean-Pierre Berenger (p.71)
    sigma_factor = (m + 1) * sigma_0 * wavelenght * omega / n_pml;
 
    // Create polynomial grading for PML regions
    for(int i = 0; i < n_pml; i++) {
        x = static_cast<double>(i) / (n_pml - 1);
        sigma = std::pow(x, m);
        
        if(direction == 'x'){
            if(upper_boundary) s = Complex(1.0, -(sigma_factor/dx[n_x - n_pml + i - 1])*sigma / (omega));
            else               s = Complex(1.0, -(sigma_factor/dx[n_pml - 1 - i])*sigma / (omega)); 

            if(n_pml > n_x) std::cout << "Error not enough points for PML" << std::endl;
            if(upper_boundary) s_profile_x[n_x - n_pml + i] = s;
            else               s_profile_x[n_pml - 1 - i] = s;
        }else if(direction == 'y'){
            if(upper_boundary) s = Complex(1.0, -(sigma_factor/dy[n_y - n_pml + i - 1])*sigma / (omega));
            else               s = Complex(1.0, -(sigma_factor/dy[n_pml - 1 - i])*sigma / (omega)); 

            if(n_pml > n_y) std::cout << "Error not enough points for PML" << std::endl;
            if(upper_boundary) s_profile_y[n_y - n_pml + i] = s;
            else               s_profile_y[n_pml - 1 - i] = s;
        }else if(direction == 'z'){
            if(upper_boundary) s = Complex(1.0, -(sigma_factor/dz[n_z - n_pml + i - 1])*sigma / (omega));
            else               s = Complex(1.0, -(sigma_factor/dz[n_pml - 1 - i])*sigma / (omega)); 

            if(n_pml > n_z) std::cout << "Error not enough points for PML" << std::endl;
            if(upper_boundary) s_profile_z[n_z - n_pml + i] = s;
            else               s_profile_z[n_pml - 1 - i] = s;                
        }else{
            std::cout << "Not one of the possible directions (x,y,z)" << std::endl;
        }
    }
}


bool CartesianCoordinateSystem::createNonUniformGrid(double refinementFactor_x, double refinementFactor_y, double refinementFactor_z, 
                                                    double x_BL, double x_BR, double x_RL,double x_RR,
                                                    double y_BL, double y_BR, double y_RL, double y_RR,
                                                    double z_BL, double z_BR, double z_RL, double z_RR,
                                                    std::vector<double> x_grid, std::vector<double> y_grid, std::vector<double> z_grid){
    
    bool non_uniform_grid = false;

    if(x_grid.size() > 0){
        non_uniform_grid = true;

        for(int i = 0; i < n_x; i++) x[i] = x_grid[i];
        for(int i=0; i<n_x-1; i++) {
            double delta = x[i + 1] - x[i];
            if (delta < 0)  throw std::runtime_error("Error: input x grid is not monotonically increasing at index " + std::to_string(i));
            dx[i] = delta;
        }
    }else if(refinementFactor_x > 0){
        non_uniform_grid = true;

        if(x_RL < x_BL || x_RR < x_RL || x_BR < x_RR) throw std::runtime_error("Error: x key grid points are not monotonically increasing: " + std::to_string(x_BL) + "," + std::to_string(x_RL) + "," + std::to_string(x_RR) + "," + std::to_string(x_BR));

        double d_normal = (length_x + x_BL - x_BR + (x_RR-x_RL)*refinementFactor_x+refinementFactor_x/(refinementFactor_x-1)*(x_RL-x_BL+x_BR-x_RR)*std::log(refinementFactor_x))/(n_x-1);
        double d_refinement = d_normal/refinementFactor_x;

        x[0] = 0.0;
        for(int i=1; i<n_x; i++){
            if(x[i-1] < x_BL){
                x[i] = x[i-1] + d_normal;
            }else if(x[i-1] < x_RL){
                double frac = (x[i-1] - x_BL)/(x_RL - x_BL);
                x[i] = x[i-1] + d_normal*(1 - frac) + d_refinement*frac;
            }else if(x[i-1] < x_RR){
                x[i] = x[i-1] + d_refinement;
            } else if(x[i-1] < x_BR){
                double frac = (x[i-1] - x_RR) / (x_BR - x_RR);
                x[i] = x[i-1] + d_refinement*(1 - frac) + d_normal*frac;
            }else{
                x[i] = x[i-1] + d_normal;
            }
        }
        double scale_domain = length_x/x[n_x-1];
        for(int i=1; i<n_x; i++) if(x[i]>x_RR) x[i] = scale_domain*x[i];

        for(int i=0; i<n_x-1; i++) dx[i] = x[i+1] - x[i];
    }

    if(y_grid.size() > 0){
        non_uniform_grid = true;

        for(int j = 0; j < n_y; j++) y[j] = y_grid[j];
        for(int j=0; j<n_y-1; j++) {
            double delta = y[j + 1] - y[j];
            if (delta < 0)  throw std::runtime_error("Error: input y grid is not monotonically increasing at index " + std::to_string(j));
            dy[j] = delta;
        }
    }else if(refinementFactor_y > 0){
        non_uniform_grid = true;

        if(y_RL < y_BL || y_RR < y_RL || y_BR < y_RR) throw std::runtime_error("Error: y key grid points are not monotonically increasing: " + std::to_string(y_BL) + "," + std::to_string(y_RL) + "," + std::to_string(y_RR) + "," + std::to_string(y_BR));

        double d_normal = (length_y + y_BL - y_BR + (y_RR-y_RL)*refinementFactor_y+refinementFactor_y/(refinementFactor_y-1)*(y_RL-y_BL+y_BR-y_RR)*std::log(refinementFactor_y))/(n_y-1);
        double d_refinement = d_normal/refinementFactor_y;

        y[0] = 0.0;
        for(int i=1; i<n_y; i++){
            if(y[i-1] < y_BL){
                y[i] = y[i-1] + d_normal;
            }else if(y[i-1] < y_RL){
                double frac = (y[i-1] - y_BL)/(y_RL - y_BL);
                y[i] = y[i-1] + d_normal*(1 - frac) + d_refinement*frac;
            }else if(y[i-1] < y_RR){
                y[i] = y[i-1] + d_refinement;
            } else if(y[i-1] < y_BR){
                double frac = (y[i-1] - y_RR) / (y_BR - y_RR);
                y[i] = y[i-1] + d_refinement*(1 - frac) + d_normal*frac;
            }else{
                y[i] = y[i-1] + d_normal;
            }
        }
        double scale_domain = length_y/y[n_y-1];
        for(int i=1; i<n_y; i++) if(y[i]>y_RR) y[i] = scale_domain*y[i];

        for(int i=0; i<n_y-1; i++) dy[i] = y[i+1] - y[i];
    }

    if(z_grid.size() > 0){
        non_uniform_grid = true;

        for(int k = 0; k < n_z; k++) z[k] = z_grid[k];
        for(int k=0; k<n_z-1; k++) {
            double delta = z[k + 1] - z[k];
            if (delta < 0)  throw std::runtime_error("Error: input z grid is not monotonically increasing at index " + std::to_string(k));
            dz[k] = delta;
        }
    }else if(refinementFactor_z > 0){
        non_uniform_grid = true;

        if(z_RL < z_BL || z_RR < z_RL || z_BR < z_RR) throw std::runtime_error("Error: z key grid points are not monotonically increasing: " + std::to_string(z_BL) + "," + std::to_string(z_RL) + "," + std::to_string(z_RR) + "," + std::to_string(z_BR));

        double d_normal = (length_z + z_BL - z_BR + (z_RR-z_RL)*refinementFactor_z+refinementFactor_z/(refinementFactor_z-1)*(z_RL-z_BL+z_BR-z_RR)*std::log(refinementFactor_z))/(n_z-1);
        double d_refinement = d_normal/refinementFactor_z;

        z[0] = 0.0;
        for(int i=1; i<n_z; i++){
            if(z[i-1] < z_BL){
                z[i] = z[i-1] + d_normal;
            }else if(z[i-1] < z_RL){
                double frac = (z[i-1] - z_BL)/(z_RL - z_BL);
                z[i] = z[i-1] + d_normal*(1 - frac) + d_refinement*frac;
            }else if(z[i-1] < z_RR){
                z[i] = z[i-1] + d_refinement;
            } else if(z[i-1] < z_BR){
                double frac = (z[i-1] - z_RR) / (z_BR - z_RR);
                z[i] = z[i-1] + d_refinement*(1 - frac) + d_normal*frac;
            }else{
                z[i] = z[i-1] + d_normal;
            }
        }
        double scale_domain = length_z/z[n_z-1];
        for(int i=1; i<n_z; i++) if(z[i]>z_RR) z[i] = scale_domain*z[i];

        for(int i=0; i<n_z-1; i++) dz[i] = z[i+1] - z[i];
    }

    return non_uniform_grid;
}


void CartesianCoordinateSystem::AddMetalPointsInVolume(Mat& matrix, Vec& b_array, std::vector<std::tuple<int, int, int>> metalLocations, double y_reflector, PetscInt low_rank, PetscInt high_rank){
    int i_idx, j_idx, k_idx;
    int i_idx_search, j_idx_search, k_idx_search;
    int pnt_idx;

    std::vector<std::vector<std::vector<bool>>> is_metal(n_x, std::vector<std::vector<bool>>(n_y, std::vector<bool>(n_z, false)));
    std::vector<std::vector<std::vector<bool>>> already_modified(n_x, std::vector<std::vector<bool>>(n_y, std::vector<bool>(n_z, false)));

    for (const auto& loc : metalLocations) is_metal[std::get<0>(loc)][std::get<1>(loc)][std::get<2>(loc)] = true;

    // Reflector plane: any point whose y-coordinate lies beyond y_reflector becomes metal too.
    for (int j = 0; j < n_y; ++j) {
        if (getY(j) > y_reflector) {
            for (int i = 0; i < n_x; ++i) {
                for (int k = 0; k < n_z; ++k) {
                    if (!is_metal[i][j][k]) {
                        is_metal[i][j][k] = true;
                        metalLocations.emplace_back(i, j, k);
                    }
                }
            }
        }
    }

    for (size_t l = 0; l < metalLocations.size(); ++l) {
        const auto& loc = metalLocations[l];

        i_idx = std::get<0>(loc);
        j_idx = std::get<1>(loc);
        k_idx = std::get<2>(loc);

        pnt_idx = point_index(i_idx, j_idx, k_idx);

        if(!(i_idx >= 0 && i_idx < n_x && j_idx >= 0 && j_idx < n_y && k_idx >= 0 && k_idx < n_z)){
            throw std::runtime_error("Metal Point is out of bounds: " + std::to_string(i_idx) + "," + std::to_string(j_idx) + "," + std::to_string(k_idx)); 
        }else{
            DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx, low_rank, high_rank);
            DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 1, low_rank, high_rank);
            DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 2, low_rank, high_rank);
        }
        
        for (int n = 0; n < 6; n++) {
            i_idx_search = i_idx; j_idx_search = j_idx; k_idx_search = k_idx;

            if(n == 0)      i_idx_search++; else if(n == 1) i_idx_search--;
            else if(n == 2) j_idx_search++; else if(n == 3) j_idx_search--;
            else if(n == 4) k_idx_search++; else if(n == 5) k_idx_search--;

            // Skip if adjacent point is out of bounds
            if (i_idx_search < 0 || i_idx_search >= n_x || j_idx_search < 0 || j_idx_search >= n_y || k_idx_search < 0 || k_idx_search >= n_z) continue;

            // Skip if adjacent is metal
            if (is_metal[i_idx_search][j_idx_search][k_idx_search]) continue;

            // Skip if already modified
            if (already_modified[i_idx_search][j_idx_search][k_idx_search]){
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx , low_rank, high_rank);
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 1, low_rank, high_rank);
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 2, low_rank, high_rank);
            }
            
            pnt_idx = point_index(i_idx_search, j_idx_search, k_idx_search);
            already_modified[i_idx_search][j_idx_search][k_idx_search] = true;

            // n = 0: +x, 1: -x, 2: +y, 3: -y, 4: +z, 5: -z
            if(n < 2){
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 1, low_rank, high_rank);
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 2, low_rank, high_rank);

                // Setting the Perpendicular component to dE_perp/dx = 0
                int row = pnt_idx;
                if(n==0 && i_idx_search+1 < n_x){ // Forward difference                    
                    if(row >= low_rank && row < high_rank){
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, -1.0, INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i_idx_search+1, j_idx_search, k_idx_search), 1.0, INSERT_VALUES);
                        VecSetValue(b_array, row, 0.0, INSERT_VALUES);
                    }
                } else if(i_idx_search > 0){ // Backward difference                    
                    if(row >= low_rank && row < high_rank){
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, 1.0, INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i_idx_search-1, j_idx_search, k_idx_search), -1.0, INSERT_VALUES);
                        VecSetValue(b_array, row, 0.0, INSERT_VALUES);
                    }
                }
            } else if(n < 4){
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx, low_rank, high_rank);
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 2, low_rank, high_rank);

                // Setting the Perpendicular component to dE_perp/dy = 0
                int row = pnt_idx + 1;
                if(n==2 && j_idx_search+1 < n_y) { // Forward difference                    
                    if(row >= low_rank && row < high_rank){
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, -1.0, INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i_idx_search, j_idx_search+1, k_idx_search) + 1, 1.0, INSERT_VALUES);
                        VecSetValue(b_array, row, 0.0, INSERT_VALUES);
                    }
                } else if(j_idx_search > 0) { // Backward difference                    
                    if(row >= low_rank && row < high_rank){
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, 1.0, INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i_idx_search, j_idx_search-1, k_idx_search) + 1, -1.0, INSERT_VALUES);
                        VecSetValue(b_array, row, 0.0, INSERT_VALUES);
                    }
                }
            } else if (n < 6) {
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx, low_rank, high_rank);
                DirichletBoundaryConditionsPoint(matrix, b_array, 0., pnt_idx + 1, low_rank, high_rank); 

                // Setting the perpendicular component to dE_perp/dz = 0
                int row = pnt_idx + 2;
                if(n==4 && k_idx_search+1 < n_z){ // Forward difference
                    if(row >= low_rank && row < high_rank){
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, -1.0, INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i_idx_search, j_idx_search, k_idx_search+1) + 2, 1.0, INSERT_VALUES);
                        VecSetValue(b_array, row, 0.0, INSERT_VALUES);
                    }
                } else if(k_idx_search > 0){ // Backward difference                    
                    if(row >= low_rank && row < high_rank){
                        SetRowZeroManually(matrix, row, n_variables, 3);
                        MatSetValue(matrix, row, row, 1.0, INSERT_VALUES);
                        MatSetValue(matrix, row, point_index(i_idx_search, j_idx_search, k_idx_search-1) + 2, -1.0, INSERT_VALUES);
                        VecSetValue(b_array, row, 0.0, INSERT_VALUES);
                    }
                }
            }          
        }          
    }
}


void CartesianCoordinateSystem::SetRowZeroManually(Mat& matrix, PetscInt row, PetscInt n_variables, int variables_per_point) {
    const int num_offsets[13][3] = {
        {0, 0, 0}, {-1, 0, 0}, {-2, 0, 0}, {1, 0, 0}, {2, 0, 0},
        {0, -1, 0}, {0, -2, 0}, {0, 1, 0}, {0, 2, 0},
        {0, 0, -1}, {0, 0, -2}, {0, 0, 1}, {0, 0, 2}
    };

    int i = getXindex(row);
    int j = getYindex(row);
    int k = getZindex(row);

    const int num_adjacent = 13 * variables_per_point;

    PetscInt* cols;
    PetscMalloc1(num_adjacent, &cols);
    PetscScalar* zeros;
    PetscMalloc1(num_adjacent, &zeros);

    PetscInt count = 0;
    for (int n = 0; n < 13; ++n) {
        int ii = i + num_offsets[n][0];
        int jj = j + num_offsets[n][1];
        int kk = k + num_offsets[n][2];

        int idx_base = point_index(ii, jj, kk);
        for (PetscInt l = 0; l < variables_per_point; ++l) {
            PetscInt col_idx = idx_base + l;
            if (col_idx < n_variables) {
                cols[count] = col_idx;
                zeros[count] = 0.0;
                ++count;
            }
        }
    }

    MatSetValues(matrix, 1, &row, count, cols, zeros, INSERT_VALUES);

    PetscErrorCode ierr;
    ierr = PetscFree(cols); CHKERRABORT(PETSC_COMM_WORLD, ierr);
    ierr = PetscFree(zeros); CHKERRABORT(PETSC_COMM_WORLD, ierr);
}


double* CartesianCoordinateSystem::normalizeDens(double* dens, int num_points, double& total_abs_power){
    double* normalized = new double[num_points];
    double total = 0;

    for (int i = 0; i < num_points; i++) total += dens[i]*getVolume(getXindex(i*3),getYindex(i*3),getZindex(i*3));

    // Printing the total density for optimization purposes
    std::ofstream fout("Outputs/p_abs.txt");
    fout << total << "\n";
    fout.close();

    total_abs_power = total;

    if (total > 0) for (int i = 0; i < num_points; ++i) normalized[i] = dens[i]/total;

    return normalized;
}


void CartesianCoordinateSystem::calculateFlux(Complex* fields, double* P_flux, double omega, int num_pts){
    // P_flux = imag(cross(E, conj(curl(E))))/(2*omega*mu_0)

    //std::vector<double> total_slice(n_y, 0.0);
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            for (int k = 0; k < n_z; k++) {
                Complex cross_x, cross_y, cross_z;
                int idx = point_index(i, j, k);
                Complex dEy_dx = 0, dEz_dx = 0, dEx_dy = 0, dEz_dy = 0, dEx_dz = 0, dEy_dz = 0;

                if(i > 0 && i < n_x-1){
                    dEy_dx = fields[point_index(i+1, j, k)+1]*(1/dx[i] - 1/(dx[i]+dx[i-1])) +
                               fields[point_index(i, j, k)+1]*(1/dx[i-1] - 1/dx[i]) +
                               fields[point_index(i-1, j, k)+1]*(1/(dx[i]+dx[i-1]) - 1/dx[i-1]);

                    dEz_dx = fields[point_index(i+1, j, k)+2]*(1/dx[i] - 1/(dx[i]+dx[i-1])) +
                               fields[point_index(i, j, k)+2]*(1/dx[i-1] - 1/dx[i]) +
                               fields[point_index(i-1, j, k)+2]*(1/(dx[i]+dx[i-1]) - 1/dx[i-1]);
                }

                if(j > 0 && j < n_y-1){
                    dEx_dy = fields[point_index(i, j+1, k)]*(1/dy[j] - 1/(dy[j]+dy[j-1])) +
                               fields[point_index(i, j, k)]*(1/dy[j-1] - 1/dy[j]) +
                               fields[point_index(i, j-1, k)]*(1/(dy[j]+dy[j-1]) - 1/dy[j-1]);
                    
                    dEz_dy = fields[point_index(i, j+1, k)+2]*(1/dy[j] - 1/(dy[j]+dy[j-1])) +
                               fields[point_index(i, j, k)+2]*(1/dy[j-1] - 1/dy[j]) +
                               fields[point_index(i, j-1, k)+2]*(1/(dy[j]+dy[j-1]) - 1/dy[j-1]);
                }

                if(k > 0 && k < n_z-1){
                    dEx_dz = fields[point_index(i, j, k+1)]*(1/dz[k] - 1/(dz[k]+dz[k-1])) +
                               fields[point_index(i, j, k)]*(1/dz[k-1] - 1/dz[k]) +
                               fields[point_index(i, j, k-1)]*(1/(dz[k]+dz[k-1]) - 1/dz[k-1]);

                    dEy_dz = fields[point_index(i, j, k+1)+1]*(1/dz[k] - 1/(dz[k]+dz[k-1])) +
                               fields[point_index(i, j, k)+1]*(1/dz[k-1] - 1/dz[k]) +
                               fields[point_index(i, j, k-1)+1]*(1/(dz[k]+dz[k-1]) - 1/dz[k-1]);
                }

                Complex curl_x_conj = std::conj(dEz_dy - dEy_dz);
                Complex curl_y_conj = std::conj(dEx_dz - dEz_dx);
                Complex curl_z_conj = std::conj(dEy_dx - dEx_dy);

                crossProduct(fields[idx], fields[idx+1], fields[idx+2], curl_x_conj, curl_y_conj, curl_z_conj, cross_x, cross_y, cross_z);

                P_flux[idx]   = std::imag(cross_x)/(2*omega*Constants::MU_0);
                P_flux[idx+1] = std::imag(cross_y)/(2*omega*Constants::MU_0);
                P_flux[idx+2] = std::imag(cross_z)/(2*omega*Constants::MU_0);
                

                //total_slice[j] += (P_flux[idx+1])*dx[i]*dz[k]; // Poynting flux in the y-direction
            }
        }
    }
    //for(int j=0; j<n_y; j++) if(j%20 == 5)std::cout << "Total Poynting Flux in slice " <<  j << ": " << total_slice[j]  << std::endl;
}


void CartesianCoordinateSystem::calculatePoynting(double* P_flux, double* absorbedPowerDensity, double* P_poynting, int num_pts){
    double total_poyting = 0;

    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            for (int k = 0; k < n_z; k++) {
                double div_P_flux_x = 0; double div_P_flux_y = 0; double div_P_flux_z = 0;
                int idx = point_index(i, j, k)/3;
                double volume = getVolume(i,j,k);

                if(i > 1 && i < n_x-1)
                    div_P_flux_x =  (P_flux[point_index(i, j, k)] - P_flux[point_index(i-1, j, k)])*cellSize(dy, j, n_y)*cellSize(dz, k, n_z);

                if(j > 1 && j < n_y-1)
                    div_P_flux_y = (P_flux[point_index(i, j, k)+1] - P_flux[point_index(i, j-1, k)+1])*cellSize(dx, i, n_x)*cellSize(dz, k, n_z);
                
                if(k > 1 && k < n_z-1)
                    div_P_flux_z = (P_flux[point_index(i, j, k)+2] - P_flux[point_index(i, j, k-1)+2])*cellSize(dx, i, n_x)*cellSize(dy, j, n_y);;

                double divergence = div_P_flux_x + div_P_flux_y + div_P_flux_z;
                
                if (n_x > 1 && (i == 0 || i == 1 || i == n_x-1)) P_poynting[idx] = 0;
                if (n_y > 1 && (j == 0 || j == 1 || j == n_y-1)) P_poynting[idx] = 0;
                if (n_z > 1 && (k == 0 || k == 1 || k == n_z-1)) P_poynting[idx] = 0;
                else P_poynting[idx] = -divergence/volume;

                //if (j < 387 or j > 665 or k < 310 or k > 570) P_poynting[idx] = 0; 
                if ( P_poynting[idx] > 6e8) P_poynting[idx] = 6e8;
                if ( P_poynting[idx] < -6e8) P_poynting[idx] = -6e8; 

                total_poyting += P_poynting[idx]*volume;
            }
        }
    }
    std::cout << "Total Poynting Power: " << total_poyting << std::endl; 

    // Printing the total density for optimization purposes
    std::ofstream fout("Outputs/power_poynting_total.txt");
    fout << total_poyting << "\n";
    fout.close();
}

void CartesianCoordinateSystem::calculatePlasmaFillingFactor(Complex* complex_conductivity, double y_center) {
    if (n_x != 1) return;

    std::vector<Complex> original_conductivity(complex_conductivity, complex_conductivity + n_y*n_z);

    const double a = length_x;  // Waveguide width in X [m]
    const int N_quad = 30;
    const double dx_quad = a / N_quad;

    // Quadrature of sin(pi*x/a) for averaging over TE_10
    const double norm_factor = M_PI / (2.0 * N_quad);   // integral of sin(pi*x/a)
    std::vector<double> x_rel_sq(N_quad);
    std::vector<double> sin_weight(N_quad);
    for (int m = 1; m < N_quad; ++m) {
        const double x_val = m * dx_quad;
        const double x_rel_m = x_val - 0.5 * a;

        x_rel_sq[m] = x_rel_m * x_rel_m;
        sin_weight[m] = std::sin(M_PI * x_val / a);
    }

    // Precomputed interpolation weights for non-uniform Y grid
    struct InterpWeight { int j0; double t; };
    struct DualInterpWeight { InterpWeight plus; InterpWeight minus; }; // Added to store both sides
    std::vector<DualInterpWeight> quad_weights(N_quad);

    // Non-uniform interpolation lookup helper
    auto get_interp_weight = [&](double y_eval) -> InterpWeight {
        if (y_eval <= y.front()) return {0, 0.0};
        if (y_eval >= y.back())  return {n_y - 2, 1.0};

        auto it = std::upper_bound(y.begin(), y.end(), y_eval);
        int j0 = std::distance(y.begin(), it) - 1;

        double y0 = y[j0];
        double y1 = y[j0 + 1];
        double t = (y_eval - y0) / (y1 - y0);

        return {j0, t};
    };

    // Main spatial iteration
    for (int j = 0; j < n_y; ++j) {
        const double y_rel = y[j] - y_center;

        if (std::abs(y_rel) > 0.014) continue; // Skipping outside tube radius

        const double y_rel_sq = y_rel * y_rel;

        // Runs once just over y to get non uniform grid weights for BOTH sides
        for (int m = 1; m < N_quad; ++m) {
            const double r_eff = std::sqrt(x_rel_sq[m] + y_rel_sq);
            
            quad_weights[m].plus = get_interp_weight(y_center + r_eff);
            quad_weights[m].minus = get_interp_weight(y_center - r_eff);
        }

        for (int k = 0; k < n_z; ++k) {
            const int idx = j * n_z + k;
            const Complex base_cond = original_conductivity[idx];

            if (std::abs(base_cond) < 1e-5) continue;

            Complex integral = 0.0;
            for (int m = 1; m < N_quad; ++m) {
                // Evaluate plus side (+r_eff)
                const int j0_p = quad_weights[m].plus.j0;
                const double t_p = quad_weights[m].plus.t;
                const Complex c0_p = original_conductivity[j0_p*n_z + k];
                const Complex c1_p = original_conductivity[(j0_p + 1)*n_z + k];
                const Complex cond_eval_p = c0_p + (c1_p - c0_p)*t_p;

                // Evaluate minus side (-r_eff)
                const int j0_m = quad_weights[m].minus.j0;
                const double t_m = quad_weights[m].minus.t;
                const Complex c0_m = original_conductivity[j0_m*n_z + k];
                const Complex c1_m = original_conductivity[(j0_m + 1)*n_z + k];
                const Complex cond_eval_m = c0_m + (c1_m - c0_m)*t_m;

                // Average the two sides together
                const Complex cond_eval_avg = 0.5 * (cond_eval_p + cond_eval_m);
                
                integral += cond_eval_avg * sin_weight[m];
            }

            const Complex averaged_cond = integral * norm_factor;
            complex_conductivity[idx] = std::sqrt(base_cond * averaged_cond);
        }
    }
}