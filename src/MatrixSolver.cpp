#include "MatrixSolver.h"

#include <iostream>
#include <fstream>

void MatrixSolver::setupGMRESSolver(Mat& systemMatrix, size_t size, PetscReal tolerance, PetscInt maxIterations, std::string preconditioner) {
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, systemMatrix, systemMatrix);
    KSPSetType(ksp, "gmres");

    PC pc;
    if(size > 1 && (preconditioner == "ilu" || preconditioner == "hypre")) preconditioner = "sor";
    KSPGetPC(ksp, &pc);
    PCSetType(pc, preconditioner.c_str());

    KSPSetFromOptions(ksp);
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(ksp, tolerance, tolerance, PETSC_DEFAULT, maxIterations);
    KSPSetUp(ksp);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}


void MatrixSolver::setupBICGSTABSolver(Mat& systemMatrix, size_t size, PetscReal tolerance, PetscInt maxIterations, std::string preconditioner) {    
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, systemMatrix, systemMatrix);
    KSPSetType(ksp, "bcgs");

    PC pc;
    if(size > 1 && (preconditioner == "ilu" || preconditioner == "hypre")) preconditioner = "bjacobi";
    KSPGetPC(ksp, &pc);
    PCSetType(pc, preconditioner.c_str());
    
    KSPSetFromOptions(ksp);
    KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    KSPSetTolerances(ksp, tolerance, tolerance, PETSC_DEFAULT, maxIterations);
    KSPSetUp(ksp);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}


void MatrixSolver::setupDirectLUSolver(Mat& systemMatrix, size_t size) {
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, systemMatrix, systemMatrix);
    KSPSetType(ksp, KSPPREONLY);
    
    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCLU);

    if (size > 1) {
        // Using SUPERLU_DIST (worse option with limited reordering)
        // PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST);
        // PetscOptionsSetValue(NULL, "-mat_superlu_dist_equil", "yes");
        // PetscOptionsSetValue(NULL, "-mat_superlu_dist_fact", "SamePattern_SameRowPerm");

        // Using MUMPS
        PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        PetscOptionsSetValue(NULL, "-mat_mumps_icntl_14", "100");
        PetscOptionsSetValue(NULL, "-mat_mumps_icntl_7", "5"); 
        PetscOptionsSetValue(NULL, "-mat_mumps_icntl_28", "2");
        PetscOptionsSetValue(NULL, "-mat_mumps_icntl_29", "2");

        // Might decrease precision
        PetscOptionsSetValue(NULL, "-mat_mumps_icntl_35", "1"); 
        PetscOptionsSetValue(NULL, "-mat_mumps_cntl_7", "5e-6"); // go to 1e-8 if solution does not look right
    }

    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
}


void MatrixSolver::setInitialGuessFromFileBinary(const std::string& filename) {
    PetscViewer viewer;
    Vec tempVec;

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer);

    VecDuplicate(solution, &tempVec);  // Ensures same layout
    VecLoad(tempVec, viewer);          // Load contents into tempVec

    VecCopy(tempVec, solution);

    VecDestroy(&tempVec);
    PetscViewerDestroy(&viewer);
}


PetscErrorCode MatrixSolver::solve(Mat& systemMatrix, Vec& rhsVector, size_t size, std::string method, std::string preconditioner) {
    if (scalar != -1) transformToScalar(scalar, systemMatrix, rhsVector);

    if (method == "GMRES") {
        sol_method = "GMRES";
        setupGMRESSolver(systemMatrix, size, tolerance, maxIterations, preconditioner);
    } else if (method == "BICGSTAB") {
        sol_method = "BICGSTAB";
        setupBICGSTABSolver(systemMatrix, size, tolerance, maxIterations, preconditioner);
    } else if (method == "direct_LU") {
        sol_method = "direct_LU";
        setupDirectLUSolver(systemMatrix, size);
    } else {
        throw std::runtime_error("Invalid solver method: " + method);
    }

    PetscErrorCode ierr = KSPSolve(ksp, rhsVector, solution);
    if (ierr) {
        PetscPrintf(PETSC_COMM_WORLD, "KSPSolve failed with error code %d\n", ierr);
    }
    
    return ierr;
}


void MatrixSolver::transformToScalar(int scalar, Mat& systemMatrix, Vec& rhsVector) {
    PetscInt rstart, rend;
    MatGetOwnershipRange(systemMatrix, &rstart, &rend);
    PetscInt local_size = rend - rstart;

    // Count how many local rows correspond to scalar-th component
    PetscInt count = 0;
    for (PetscInt i = 0; i < local_size; ++i) {
        PetscInt global_idx = rstart + i;
        if (global_idx % 3 == scalar) count++;
    }

    PetscInt* scalar_indices;
    PetscMalloc1(count, &scalar_indices);

    // Fill scalar_indices with global indices matching scalar-th component
    PetscInt idx = 0;
    for (PetscInt i = 0; i < local_size; ++i) {
        PetscInt global_idx = rstart + i;
        if (global_idx % 3 == scalar) {
            scalar_indices[idx++] = global_idx;
        }
    }

    IS is_scalar;
    ISCreateGeneral(PETSC_COMM_WORLD, count, scalar_indices, PETSC_COPY_VALUES, &is_scalar);
    //PetscFree(scalar_indices);

    Mat systemMaxwell_scalar;
    MatCreateSubMatrix(systemMatrix, is_scalar, is_scalar, MAT_INITIAL_MATRIX, &systemMaxwell_scalar);

    Vec b_vector_scalar;
    VecGetSubVector(rhsVector, is_scalar, &b_vector_scalar);

    Mat oldMat = systemMatrix;
    Vec oldVec = rhsVector;

    systemMatrix = systemMaxwell_scalar;

    // create a true copy of the subvector
    VecDuplicate(b_vector_scalar, &rhsVector);
    VecCopy(b_vector_scalar, rhsVector);

    VecRestoreSubVector(oldVec, is_scalar, &b_vector_scalar);

    MatDestroy(&oldMat);
    VecDestroy(&oldVec);
    ISDestroy(&is_scalar);
}


void MatrixSolver::printConvergenceInfo() {
    if(sol_method != "direct_LU"){
        PetscInt num_iterations;
        KSPGetIterationNumber(ksp, &num_iterations);
        PetscPrintf(PETSC_COMM_WORLD, "Solver iterations: %d\n", num_iterations);

        PetscReal residual_norm;
        KSPGetResidualNorm(ksp, &residual_norm);
        PetscPrintf(PETSC_COMM_WORLD, "Final Residual Norm: %g\n", residual_norm);
    }
}