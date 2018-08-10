

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "../std_math/sparse_matrices.hpp"
#include "../std_math/matrix_generators.hpp"
#include "../std_math/std_math.hpp"

using namespace std_math;
// pointer to dynamic storage scheme
// =================================
dynamic_storage_scheme<double> *DSS;

int main()
{

    {
    // creation of empty input matrix
    // ==============================
    input_storage_scheme<double> iss;
    load_scheme_from_file("input_matrices/input_matrix_sor2.txt", &iss);
    input_storage_scheme<double> ISS(iss);

    iss = ISS;
    // dimension of the linear problem
    // ===============================
    const int Nsize = ISS.order;

    try
    {
        DSS = new dynamic_storage_scheme<double>(ISS, 1, 1);
    }
    catch (std::bad_alloc error)
    {
        std::cout << error.what() << std::endl;
        system("PAUSE");
        return 0;
    }
    DSS->print_scheme_to_file("output_matrices/out_matrix_1.mxout");

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

    std_math::load_vector(b, Nsize, "saved_vectors/b2.vtr");

    for (int i = 0; i < Nsize; i++)
        x[i] = 0.0f;

    //DSS->print_scheme_to_file("output_matrices/out_matrix_1.mxout");
    // prepare matrix to SOR itarations
    // ================================
    DSS->print_sparsity_pattern("output_matrices/out_pattern_1.patt");
    DSS->iterative_preparation();
    //DSS->print_scheme_to_file("output_matrices/out_matrix_2.mxout");


    DSS->iterative_refinement(ISS, SOR, x, b, 0.0000000000001, 100, d, r, y);


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
