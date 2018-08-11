#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "../../sparse_matrices.hpp"

int main()
{
    // declaration of the using namespace
    using namespace std_math;

    // declaration of two input schemes
    input_storage_scheme<double> ISS1;

    // first way of loading scheme from file
    load_scheme_from_file("example_scheme.txt", &ISS1);

    // second way of loading scheme from file
    input_storage_scheme<double> ISS2(ISS1);

    ISS2.add_element(1,2,3);

    // print schemes to console
    std::cout << ISS1;
    std::cout << ISS2;

    system("pause");

    return 0;



}