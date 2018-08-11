#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <complex>
#include "../../sparse_matrices.hpp"

// used namespaces
// ===============
using namespace std;

// function used for creation random complex numbers
// =================================================
complex<double> get_complex_random_number()
{
const int freq = 2;

    int rand_val1 = rand();
    int rand_val2 = rand();
    if (rand_val1 % freq == 0)
    {
        complex<double> val(static_cast<double>(rand_val1 % 10) ,
                            static_cast<double>(rand_val2 % 10));
        return val;
    }
    else
        return complex<double>(0,0);
}

const int N = 5;
int main()
{
    // initialize random seed
    // ======================
    srand (time(NULL));

    // vectors of the equation Ax=b
    // ============================
    complex<double> x[N], b[N];

    // declaration of the input scheme
    // ===============================
    input_storage_scheme<complex<double> > ISS(N,N);

    // initialize of the input scheme and vector b
    // ===========================================
    for (int row = 0; row < N; ++row)
    {
        for (int col = 0; col < N; ++col)
        {
            const complex<double> val = get_complex_random_number();
            if (val.real() != 0 || val.imag() != 0)
                ISS.add_element(val, row, col);
        }
        b[row] = get_complex_random_number();
    }

    // printing vector b
    // =================
    cout << endl << "b:\n";
    for (int i = 0; i < N; ++i) cout << b[i] << endl;

    // creation of the dynami scheme
    // =============================
    dynamic_storage_scheme<complex<double> > DSS(ISS, 3, 0.7, ROL_INIT);

    // printing dynamic scheme in its initial form
    // ===========================================
    cout << "\ndynamic_storage_scheme in its initial form\n";
    cout << DSS;

    // attempt to decompose the matrix A
    // =================================
    try
    {
        DSS.LU_decomposition(MARKOWITZ_COST,10,2,0.001,true);
    }
    catch (std::exception e)
    {
        std::cout << e.what();
        return 0;
    }

    // printing dynamic scheme after Gauss elimination
    // ===============================================
    cout << "\n\ndynamic_storage_scheme after GE decomposition\n";
    cout << DSS;

    // obtaining inital solution
    // =========================
    DSS.solve_LU(x,b);

    // performing interative refinement
    // ================================
    DSS.iterative_refinement(ISS,x,b,0.00000000001,10);

    // printing solution vector x
    // ==========================
    cout << "\nx:\n";
    for (int i = 0; i < N; ++i)
        cout << x[i] << endl;

    // this function multiply obtained solution x by input matrix and sore it to vector b
    // ==================================================================================
    CheckAxb(ISS, x, b);

    // printing check vector b
    // =======================
    cout << "\ncheck b:\n";
    for (int i = 0; i < N; ++i)
        cout << b[i] << endl;

    // halt program
    // ============
    system("pause");
}
