#ifndef MATRIX_GENERATORS
#define MATRIX_GENERATORS

#include <stdlib.h>
#include <time.h>
#include <stdexcept>

class matrix_generator
{
private:
    // sizes of generated matrix
    size_t number_of_rows, number_of_columns;
    // higher of the above values
    size_t N;
    // table of indexes from 0 to max(number_of_rows, number_of_columns) for performing non singular matrix
    size_t *permutation_table;
    // disable default constructor
    matrix_generator(void);

public:
    // constructor of matrix_generator
    matrix_generator(size_t n_of_rows, size_t n_of_cols);
    // destructor
    ~matrix_generator(void)
    {
        delete[] permutation_table;
    }
    // randomize permutation table
    void randomize_permutation_table(void);
    // method sets permutation table do identity
    void set_permutation_table_to_identity(void);
    // method returns random element
    double get_random_double_in_range(size_t range, size_t freq, size_t row, size_t col);
};





#endif