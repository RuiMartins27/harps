#ifndef SOLVER_H  // Include guard to prevent double inclusion
#define SOLVER_H

#include <petscksp.h>

class MatrixSolver {
private:
    KSP ksp;
    Vec solution;

    PetscInt num_variables;
    PetscReal tolerance;
    PetscInt maxIterations;
    std::string sol_method;

    int scalar;  // -1 for vector system, 0 for Ex, 1 for Ey, 2 for Ez

    void setupGMRESSolver(Mat& systemMatrix, size_t size, PetscReal tolerance, PetscInt maxIterations, std::string preconditioner = "ilu");
    void setupBICGSTABSolver(Mat& systemMatrix, size_t size, PetscReal tolerance, PetscInt maxIterations, std::string preconditioner = "ilu");
    void setupDirectLUSolver(Mat& systemMatrix, size_t size);

public:
    MatrixSolver(PetscInt numVariables, PetscReal tolerance, PetscInt maxIterations, int scalar = -1)
                        : num_variables(numVariables), tolerance(tolerance), maxIterations(maxIterations), scalar(scalar) {
        
        if(scalar != -1 && scalar != 0 && scalar != 1 && scalar != 2) throw std::invalid_argument("Scalar must be -1 (vector system), 0 (Ex), 1 (Ey), or 2 (Ez)");
        
        if(scalar != -1) num_variables /= 3;

        VecCreate(PETSC_COMM_WORLD, &solution);
        VecSetSizes(solution, PETSC_DECIDE, num_variables);
        VecSetFromOptions(solution);
        
        if(scalar == -1) VecSetBlockSize(solution, 3);
        
        //PetscOptionsSetValue(NULL, "-ksp_monitor_true_residual", "");
        //PetscOptionsSetValue(NULL, "-ksp_converged_reason", "");
        //PetscOptionsSetValue(NULL, "-ksp_view", "");
    }

    ~MatrixSolver() {
        VecDestroy(&solution);
        KSPDestroy(&ksp);
    }
    
    void setInitialGuessFromFileBinary(const std::string& filename);
    void transformToScalar(int scalar, Mat& systemMatrix, Vec& rhsVector);
    PetscErrorCode solve(Mat& systemMatrix, Vec& rhsVector, size_t size, std::string method = "direct_LU", std::string preconditioner = "ilu");
    Vec getSolution() const { return solution; }
    
    void printConvergenceInfo();
};

#endif // SOLVER_H
