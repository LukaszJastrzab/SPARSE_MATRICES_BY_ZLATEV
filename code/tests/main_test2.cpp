

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "std_math\matrix_storage_schemes.hpp"
#include "std_math\matrix_generators.hpp"

// dimension of the linear problem
// ===============================
const int Nsize = 150;

// pointer to dynamic storage scheme
// =================================
dynamic_storage_scheme<double> *DSS;

int main()
{

    {
    // create a random matrix generator to generate lerger random non-singular matrix
    // ==============================================================================
    matrix_generator gen(Nsize, Nsize);

    // initiate random number generator
    // ================================
    srand((size_t)time(NULL));

    // creation of empty input matrix
    // ==============================
    input_storage_scheme<double> ISS(Nsize, Nsize);

    // create random test matrix
    // =========================
    matrix_generator m_gen(Nsize, Nsize);
    for (size_t i = 0; i < Nsize; ++i)
        for (size_t j = 0; j < Nsize; ++j)
        {
            double val = m_gen.get_random_double_in_range(20, 2, i, j);
            if (val != 0)
                ISS.add_element(val, i, j);
        }

    permute_input_matrix_elements_test(&ISS);

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
        DSS = new dynamic_storage_scheme<double>(ISS, 500, 0.8);
    }
    catch (std::bad_alloc error)
    {
        std::cout << error.what() << std::endl;
        system("PAUSE");
        return 0;
    }

    // print sparsity patern of generated matrix to file
    // =================================================
    DSS->print_sparsity_pattern("output_matrices/out_pattern_1.mxout");

    // prepare matrix to SOR itarations
    // ================================
    DSS->GAUSS_decomposition(MARKOWITZ_COST, 30, 50, 0.000000001);
    DSS->solve_LU(x, b, y);
    DSS->iterative_refinement(ISS, LU, x, b, 0.0000000000001, 10, d, r, y);

    // print sparsity patern of decomposed matrix to file
    // ==================================================
    DSS->print_sparsity_pattern("output_matrices/out_pattern_2.mxout");
    system("PAUSE");


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
    ISS.count_rasidual_vector(x,b,r);
    for (int i = 0; i < Nsize; i++)
        std::cout << std::fixed << std::setw(12) << r[i] << " ";

    std::cout << std::endl;
    system("PAUSE");

    // delete allocated memory
    // =======================
    delete DSS;
    delete[] x;
    delete[] b;
    delete[] r;
    delete[] d;
    delete[] y;
    }


system("PAUSE");
return 0;
}
