

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "MatrixStorageSchemes.hpp"

// dimension of the linear problem
// ===============================
const int Nsize = 8;

// pointer to dynamic storage scheme
// =================================
dynamic_storage_scheme<double> *M;

int main()
{

    {
    system("PAUSE");

    // initiate random number generator
    // ================================
    srand((size_t)time(NULL));

    // creation of empty input matrix
    // ==============================
    input_storage_scheme<double> Min(Nsize, Nsize);

    // loading scheme from the file
    // ============================
    load_scheme_from_file("input_error.txt", Min);

    // print stored matrix
    // ===================
    std::cout << Min;

    // allocation of vectors
    // =====================
    double *b = new double[Nsize];  // vector b from Ax = b
    double *x = new double[Nsize];  // searching solution x = A^{-1}b

    // vectors below are only optional, created for savings in process performance time
    // ================================================================================
    double *r = new double[Nsize]; // memory allocated for residual vector in performance of iterative refinement
    double *d = new double[Nsize]; // memory allocated for added vector in performance of iterative refinement
    double *y = new double[Nsize]; // memory allocated for helpfull vector in solving LUx = b (Ly = b, Ux = y)

    // fill vector b with random numbers
    // =================================
    for (int i = 0; i < Nsize; i++)
        b[i] = static_cast<double>(rand() % 1000 + 1);

    // creation of dynamic storage scheme including exception handling
    // ===============================================================
    try
    {
        M = new dynamic_storage_scheme<double>(Min, 200.5, 0.8);
    }
    catch (std::bad_alloc error)
    {
        std::cout << error.what() << std::endl;
        system("PAUSE");
        return 0;
    }

    // print dynamic storage scheme to file
    // ====================================
    M->print_scheme_to_file("main_matrix_1.mxout");

    // print sparsity pattern to file
    // ==============================
    M->print_sparsity_pattern("sparsity_patern_before.pat");

    // decomposition of a matrix including exception handling
    // ======================================================
    try
    {
        M->GAUSS_decomposition(MARKOWITZ_COST, 8, 5, 0.000000001, true);
    }
    catch (std::exception e)
    {
        std::cout << e.what() << std::endl << M->get_LogErrors();
        system("PAUSE");

        delete M;
        delete[] x;
        delete[] b;
        delete[] r;
        delete[] d;
        delete[] y;
        return 0;
    }

    system("PAUSE");

    // print decomposed matrix to file
    // ===============================
    M->print_scheme_to_file("main_matrix_2.mxout");

    // print sparsity pattern of decomposed matrix to file
    // ===================================================
    M->print_sparsity_pattern("sparsity_patern_after.pat");

    // solve the equation to obtain first aproximation of the solution
    // ===============================================================
    M->solve_LU(x, b);

    // improve first aproximation
    // ==========================
    M->iterative_refinement(Min, x, b, 0.0000000000000000001, 20, d, r);

    // print right side vector b
    // =========================
    std::cout << std::endl << std::endl;
    std::cout << "b:" << std::endl;
    for (int i = 0; i < Nsize; i++)
        std::cout << std::fixed << std::setw(12) << b[i] << " ";

    // print counted solution
    // ======================
    std::cout << "\nx:" << std::endl;
    for (int i = 0; i < Nsize; i++)
        std::cout << std::fixed << std::setw(12) << x[i] << " ";

    // print residual vector of the counted solution
    // ==============================================
    std::cout << "\nr:" << std::endl;
    Min.count_rasidual_vector(x,b,r);
    for (int i = 0; i < Nsize; i++)
        std::cout << std::fixed << std::setw(12) << r[i] << " ";

    std::cout << std::endl;
    system("PAUSE");

    // delete allocated memory
    // =======================
    delete M;
    delete[] x;
    delete[] b;
    delete[] r;
    delete[] d;
    delete[] y;
    }


system("PAUSE");
return 0;
}
