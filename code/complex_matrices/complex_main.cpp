#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "../std_math/sparse_matrices.hpp"
#include "../std_math/std_math.hpp"

using namespace std;
using namespace std_math;

complex_number<double> get_random_element(void)
{
    complex_number<double> val(static_cast<double>(rand()) / 100, static_cast<double>(rand()) / 100);
    return val;
}

int main()
{
    srand((size_t)time(NULL));
    const size_t N = 5;
    input_storage_scheme<complex_number<double> > ISS(N,N);

    for (int row = 0; row < N; ++row)
        for (int col = 0; col < N; ++col)
        {
            if (rand() % 2 == 0)
                ISS.add_element(get_random_element(),row,col);
        }

    dynamic_storage_scheme<complex_number<double> > DSS(ISS,5,0.7,ROL_INIT);
    try
    {
        DSS.LU_decomposition(MARKOWITZ_COST,10,2,0.001,true);
    }
    catch (std::exception e)
    {
        std::cout << e.what();
    }

    return 0;
}