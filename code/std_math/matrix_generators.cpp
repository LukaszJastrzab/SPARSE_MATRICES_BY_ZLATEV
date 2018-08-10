#include "matrix_generators.hpp"


//----------------------------------------------------------------------------------------- matrix_generator
/**
*  Construtor of matrix generator
*
*  @param n_of_rows         - number of rows of a matrix we want to generate
*  @param n_of_cols         - number of columns of a matrix we want to generate
*/
//----------------------------------------------------------------------------------------------------------
matrix_generator::matrix_generator(size_t n_of_rows, size_t n_of_cols)
 :
 number_of_rows(n_of_rows),
 number_of_columns(n_of_cols),
 N(n_of_rows > n_of_cols ? n_of_rows : n_of_cols),
 permutation_table(NULL)
{
    permutation_table = new size_t[N];
    for (size_t i = 0; i < N; ++i)
        permutation_table[i] = i;
}

//------------------------------------------------------------------------------ randomize_permutation_table
/**
*  Function randomly changing elements of a permutation_table
*/
//----------------------------------------------------------------------------------------------------------
void matrix_generator::randomize_permutation_table(void)
{
    srand((size_t)time(NULL));
    for (size_t i = 0; i < 5 * N; ++i)
    {
        size_t idx1 = rand() % N;
        size_t idx2 = rand() % N;
        if (idx1 != idx2)
        {
            size_t val = permutation_table[idx1];
            permutation_table[idx1] = permutation_table[idx2];
            permutation_table[idx2] = val;
        }
    }
}
//------------------------------------------------------------------------------ randomize_permutation_table
/**
*  Method sets permutation teble to identity
*/
//----------------------------------------------------------------------------------------------------------
void matrix_generator::set_permutation_table_to_identity(void)
{
    for (size_t i = 0; i < N; ++i)
        permutation_table[i] = i;
}

//---------------------------------------------------------------------------------- random_element_in_range
/**
*  Method is used to create a random matrix which each of its row and each column has at least one 
*  nonzero element, this ensure non-singular patrix in ~99% other non-zeroes can lie only in
*  specified range form diagonal with specified frequently
*
*  @param range        - specified range from diagonal where non-zeros can appear
*  @param freq         - specified the frequency of nonzero elements
*/
//----------------------------------------------------------------------------------------------------------
double matrix_generator::get_random_double_in_range(size_t range, size_t freq, size_t row, size_t col)
{
    if (row == 0 && col == 0)
        srand((size_t)time(NULL));

    if (row > number_of_rows || col > number_of_columns)
        throw std::out_of_range("get_random_double_in_range: out_of_range");

    if (permutation_table[row] == col || (abs((int)row - (int)col) <= (int)range && rand() % freq == 0))
        return static_cast<double>(rand()/100.0f);

return 0.0f;
}


