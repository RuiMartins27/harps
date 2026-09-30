#include <petscksp.h>
#include <cmath>
#include <iostream>

int main(int argc, char* argv[]) {
    
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    int rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    
    int num_variables = 18;
    Vec b_vector;
    
    VecCreate(PETSC_COMM_WORLD, &b_vector);
    VecSetSizes(b_vector, PETSC_DECIDE, num_variables);
    VecSetFromOptions(b_vector);
    
    // Initialize all values to 0.0
    VecSet(b_vector, 0.01);
    
    // Get the local range for this process
    PetscInt low, high;
    VecGetOwnershipRange(b_vector, &low, &high);
    std::cout << "Rank " << rank << " owns rows " << low << " to " << high << std::endl;
    
    for (int i = 0; i < num_variables; i++) {
        if (i % 1 == 0) {
            PetscReal value = i+2;
            if(i >= low && i < high) VecSetValue(b_vector, i, value, INSERT_VALUES);
        }
    }
    for (int i = 0; i < num_variables; i++) {
        if (i % 2 == 0) {
            PetscReal value = -i * 10;
            if(i >= low && i < high) VecSetValue(b_vector, i, value, INSERT_VALUES);
        }
    }
    for (int i = 0; i < num_variables; i++) {
        if (i % 2 == 0) {
            PetscReal value = -i * 0.1;
            if(i >= low && i < high) VecSetValue(b_vector, i, value, INSERT_VALUES);
        }
    }
    for (int i = 0; i < num_variables; i++) {
        if (i % 2 == 0) {
            PetscReal value = -i * 12345687;
            if(i >= low && i < high) VecSetValue(b_vector, i, value, INSERT_VALUES);
        }
    }
    // Assemble the vector
    VecAssemblyBegin(b_vector);
    VecAssemblyEnd(b_vector);

    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, num_variables, num_variables);
    MatSetFromOptions(A);
    MatSetOption(A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
    MatSetUp(A);

    for(int i = 0; i < num_variables; i++) {
        for (int j = 0; j < num_variables; j++) {
            if (i == j) {
                PetscReal value = 2.0;
                MatSetValue(A, i, j, value, INSERT_VALUES);
            }
            if(i > 0) {
                PetscReal value = -1.0;
                if (j== i-1) MatSetValue(A, i, j, value, INSERT_VALUES);
            }
            if(i < num_variables - 1) {
                PetscReal value = -1.0;
                if (j== i+1) MatSetValue(A, i, j, value, INSERT_VALUES);
            }
        }
    }

    MatSetValue(A, 2, 10, 27, INSERT_VALUES);


    PetscInt row;

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    row = num_variables/2;
    MatZeroRows(A, 1, &row, 1.0, PETSC_NULLPTR, PETSC_NULLPTR);


    MatSetValue(A, 3, 11, 17, INSERT_VALUES);


    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    row = row + 1;
    MatZeroRows(A, 1, &row, 1.0, PETSC_NULLPTR, PETSC_NULLPTR);

    MatSetValue(A, 4, 12, 7, INSERT_VALUES);

    

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    PetscInt       *indices;
    PetscErrorCode ierr;
    ierr = PetscMalloc1(1, &indices); CHKERRQ(ierr);

    PetscInt row_zero = 5;
    indices[0] = row_zero;

    //ierr = MatZeroRowsColumns(A, 1, indices, 1.0, PETSC_NULLPTR, PETSC_NULLPTR); CHKERRQ(ierr);
    ierr = MatZeroRows(A, 1, &row_zero, 1., PETSC_NULLPTR, PETSC_NULLPTR);

    ierr = PetscFree(indices); CHKERRQ(ierr);




    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    // View the vector and the matrix
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "vector_test", &viewer);
    VecView(b_vector, viewer);
    PetscViewerDestroy(&viewer);

    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "matrix_test", &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);
    
    PetscFinalize();
    return 0;
}