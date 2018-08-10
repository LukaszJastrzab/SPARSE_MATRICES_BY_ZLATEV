

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "../std_math/sparse_matrices.hpp"



double get_random_element(void)
{
    int rand_val = rand();
    if (rand_val % 10 == 0)
        return static_cast<double>(rand_val) / 100;
    else
        return 0;
}


int main()
{

    srand((size_t)time(NULL));             // initialization of randomizer
    const size_t N = 100;                  // size of the problem
    double x[N], b[N];                     // vectors of the equation Ax=b
    std_math::input_storage_scheme<double> ISS(N,N); // input scheme containing matrix A

    // random initialization of the matrix A and vektor b
    for (int row = 0; row < N; ++row)
    {
        for (int col = 0; col < N; ++col)
        {
            double val = get_random_element();
            if (val != 0)
                ISS.add_element(val,row,col);
        }
        b[row] = static_cast<double>(rand()) / 100;
    }

    // creation of the dynamic scheme containing matrix A
    std_math::dynamic_storage_scheme<double> DSS(ISS,5,0.7,std_math::ROL_INIT);

    // attempt to decompose the matrix A
    try
    {
        DSS.LU_decomposition(std_math::MARKOWITZ_COST,10,2,0.001,true);
    }
    catch (std::exception e)
    {
        std::cout << e.what();
        return 0;
    }

    // solution of linear system with iterative improvement
    DSS.solve_LU(x,b);
    DSS.iterative_refinement(ISS,x,b,0.00000000001,10);

return 0;
}
