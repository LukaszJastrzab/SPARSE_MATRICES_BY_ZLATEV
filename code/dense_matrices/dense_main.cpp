
#include <iostream>
#include "../std_math/dense_matrices.hpp"
#include "../std_math/std_math.hpp"

using namespace std_math;

int main()
{
    dense_matrix<double> A;
    load_matrix_from_file("input_matrices/input_matrix_sor2.txt", &A);
    std::cout << A << std::endl;

    const size_t Nsize = A.order;
    double *b = new double[Nsize];
    std_math::load_vector(b, Nsize, "saved_vectors/b2.vtr");
    double *x = new double[Nsize];
    for (size_t i = 0; i < Nsize; i++)
        x[i] = 0.0f;

    A.iterative_refinement(SOR, x, b, 0.0000000001, 100);

    system("pause");
    delete[] b;
    delete[] x;



    return 0;
}