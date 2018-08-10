#ifndef DENSE_MATRICES
#define DENSE_MATRICES

#include <stdlib.h>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <fstream>
#include <iomanip>

#include <std_math\std_math.hpp>

namespace std_math{

//-----------------------------------------------------------------------------------------------//
//                                    dense_matrx class                                          //
//-----------------------------------------------------------------------------------------------//
template <typename TYPE>
class dense_matrix
{
public:
    /// basic data which describes the matrix
    const size_t number_of_rows,
                 number_of_columns,
                 order;
    /// array of elements
    TYPE** a;

    /// relaxation parameter
    double omega;
public:
    /// default and convertion constructor
    dense_matrix(TYPE val = 0);
    /// main constructor
    dense_matrix(size_t n_of_rows, size_t n_of_cols, bool fill = false, TYPE val = 0);
    /// destructor
    ~dense_matrix();

    /// method puts the value in specified place in matrix
    void set_value(TYPE val, size_t row, size_t col)
    {
        if (row >= number_of_rows || col >= number_of_columns)
            throw std::out_of_range("dense_matrix<TYPE>::set_value: out_of_range");
        else
            a[row][col] = val;
    }
    /// method clear everything or sets specific matrix
    void clear(bool all = true, TYPE val = TYPE(0), bool diag = false);
    /// methode resizes the matrix without saving data
    void resize(size_t n_of_rows, size_t n_of_cols);

    /// method counts residual vector r_i = Ax_i - b
    void count_rasidual_vector(TYPE *x, TYPE *b, TYPE *r) const;
    /// method performs one iteration of the Jackobi method
    void J_iteration(TYPE *x = NULL, TYPE *b = NULL, TYPE *prev_x = NULL) const;
    /// method performs one iteration of the Gauss-Seidel method
    void GS_iteration(TYPE *x = NULL, TYPE *b = NULL, TYPE *prev_x = NULL) const;
    /// method performs one iteration of SOR method
    void SOR_iteration(TYPE *x, TYPE *b, TYPE *prev_x) const;
    /// Method improves the accuracy of the solution
    void iterative_refinement(const EQUATION_METHOD method, TYPE *x, TYPE *b, double acc, size_t max_it, const dense_matrix<TYPE>* DS = NULL, TYPE *d = NULL, TYPE *r = NULL, TYPE *y = NULL)const;

private:
    /// method for loading matrix from 
    template <typename TYPE2>
    friend void load_matrix_from_file(const char* file_name, dense_matrix<TYPE2>* DM);
    /// standard output operator
    template <typename TYPE2>
    friend std::ostream& operator<< (std::ostream& out, const dense_matrix<TYPE2>& DM);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//                                  DEFINITIONS OF METHODS                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------------------------------------- clear
/// statndard constructor
///
/// @param all        - flag indicating if memory should be relesed
/// @param val        - if memory shouldn't be relesed then fill it with this value
/// @param diag       - fill with val all or only diagonal (zero the rest)
///
//--------------------------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::clear(bool all, TYPE val, bool diag)
{
    if (all)
    {
        for (size_t i = 0; i < number_of_rows; ++i)
            delete[] a[i];
        delete[] a;
        number_of_rows = number_of_columns = 0;
    }
    else
    {
        for (size_t i = 0; i < number_of_rows; ++i)
            for (size_t j = 0; j < number_of_columns; ++j)
            {
                if ((i != j && !diag ) || i == j)
                    a[i][j] = val;
                else
                    a[i][j] = TYPE(0);
            }
    }
}
//------------------------------------------------------------------------------------------ resize
/// methods resizes the array (don't keep previous values)
///
/// @param n_of_rows        - number of rows of new size
/// @param n_of_cols        - number of columns of new size
///
/// @throw bad_alloc        - no enough memory
//--------------------------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::resize(size_t n_of_rows, size_t n_of_cols)
{
    if (n_of_rows != number_of_rows || n_of_cols != number_of_columns)
    {
        for (size_t i = 0; i < number_of_rows; ++i)
            delete[] a[i];
        delete[] a;
        a = NULL;

        const_cast<size_t&>(number_of_rows) = n_of_rows;
        const_cast<size_t&>(number_of_columns) = n_of_cols;
        const_cast<size_t&>(order) = (n_of_cols < n_of_rows ? n_of_cols : n_of_rows);

        try
        {
            a = new TYPE*[n_of_rows];
            for (size_t i = 0; i < n_of_rows; ++i)
                a[i] = NULL;
            for (size_t i = 0; i < n_of_rows; ++i)
                a[i] = new TYPE[n_of_cols];
        }
        catch (std::bad_alloc)
        {
            if (a != NULL)
                for (size_t i = 0; i < n_of_rows; ++i)
                    delete[] a[i];
            delete[] a;
            throw std::bad_alloc("dense_matrix<TYPE>::resize: bad_alloc");
        }
    }
}

//------------------------------------------------------------ count_rasidual_vector
/**
*  Method used to counting residual vector
*
*  @param x                - [in]  obtained solution
*  @param b                - [in]  right sight vestor of equation Ax=b
*  @param r                - [out] counting residual vestor
*/
//----------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::count_rasidual_vector( TYPE *x,
                                                TYPE *b,
                                                TYPE *r
                                                ) const
{
    for (size_t row = 0; row < number_of_rows; ++row)
    {
        r[row] = -b[row];
        for (size_t col = 0; col < number_of_columns; ++col)
            r[row] += (x[col] * a[row][col]);
    }
}

//------------------------------------------------------------------------------------ dense_matrix
/// statndard constructor
///
/// @param n_of_rows        - number of rows of stored matrix
/// @param n_of_cols        - number of columns of stored matrix
///
/// @throw bad_alloc
//--------------------------------------------------------------------------------------------------
template <typename TYPE>
dense_matrix<TYPE>::dense_matrix(TYPE value)
:
number_of_rows(1),
number_of_columns(1),
order(1),
a(NULL),
omega(1)
{
    a = new TYPE*[1];
    a[0] = new TYPE[1];
    a[0][0] = value;
}

//------------------------------------------------------------------------------------ dense_matrix
/// statndard constructor
///
/// @param n_of_rows        - number of rows of stored matrix
/// @param n_of_cols        - number of columns of stored matrix
/// @param fill             - fillin flag (default false)
/// @param val              - fillin value
///
/// @throw bad_alloc
//--------------------------------------------------------------------------------------------------
template <typename TYPE>
dense_matrix<TYPE>::dense_matrix( size_t n_of_rows,
                                  size_t n_of_cols,
                                  bool fill,
                                  TYPE val
                                  )
:
number_of_rows(n_of_rows),
number_of_columns(n_of_cols),
order(n_of_rows < n_of_cols ? n_of_rows : n_of_cols),
a(NULL),
omega(1)
{
    try
    {
        a = new TYPE* [n_of_rows];
        for (size_t i = 0; i < n_of_rows; ++i)
            a[i] = NULL;
        for (size_t i = 0; i < n_of_rows; ++i)
            a[i] = new TYPE[n_of_cols];
    }
    catch (std::bad_alloc)
    {
        if (a != NULL)
            for (size_t i = 0; i < n_of_rows; ++i)
                delete[] a[i];
        delete[] a;
        throw std::bad_alloc("dense_matrix: dense_matrix bad_alloc ");
    }
    if (fill)
    {
        for (size_t i = 0; i < n_of_rows; ++i)
            for (size_t j = 0; j < n_of_cols; ++j)
                a[i][j] = val;
    }
}
//----------------------------------------------------------------------------------- ~dense_matrix
/// destructor
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
dense_matrix<TYPE>::~dense_matrix()
{
    for (size_t i = 0; i < number_of_rows; ++i)
        delete[] a[i];
    delete[] a;
}
//------------------------------------------------------------------------------------- J_iteration
/// method performs one iteration of Jacobi method
///
/// @param x            - solution
/// @param b            - right side vector b of equation Ax=b
/// @param prev_x       - previous appoximation (or starting point)
///
/// @throw excpetion    - varius reasons
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::J_iteration( TYPE *x,
                                      TYPE *b,
                                      TYPE *prev_x
                                      ) const
{
    if (number_of_rows != number_of_columns)
        throw std::exception("dense_matrix<TYPE>::J_iteration: matrix is not squered");
    if (!x || !b || !prev_x)
        throw std::exception("dense_matrix<TYPE>::J_iteration: invalid parameter");

    for (size_t row = 0; row < order; ++row)
    {
        x[row] = b[row];
        for (size_t col = 0; col < order; ++col)
            if (row != col)
                x[row] -= (a[row][col] * prev_x[col]);

        x[row] /= a[row][row];
    }
}

//------------------------------------------------------------------------------------ GS_iteration
/// method performs one iteration of Gauss-Seidel method
///
/// @param x            - solution
/// @param b            - right side vector b of equation Ax=b
/// @param prev_x       - previous appoximation (or starting point)
///
/// @throw excpetion    - varius reasons
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::GS_iteration( TYPE *x,
                                       TYPE *b,
                                       TYPE *prev_x
                                       ) const
{
    if (number_of_rows != number_of_columns)
        throw std::exception("dense_matrix<TYPE>::GS_iteration: matrix is not squered");
    if (!x || !b || !prev_x)
        throw std::exception("dense_matrix<TYPE>::GS_iteration: invalid parameter");

    for (size_t row = 0; row < order; ++row)
    {
        x[row] = b[row];
        for (size_t col = 0; col < row; ++col)
            x[row] -= (a[row][col] * x[col]);
        for (size_t col = row + 1; col < order; ++col)
            x[row] -= (a[row][col] * prev_x[col]);
        x[row] /= a[row][row];
    }
}
//------------------------------------------------------------------------------------ GS_iteration
/// method performs one iteration of the SOR method
///
/// @param x            - solution
/// @param b            - right side vector b of equation Ax=b
/// @param prev_x       - previous appoximation (or starting point)
///
/// @throw excpetion    - varius reasons
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::SOR_iteration( TYPE *x,
                                        TYPE *b,
                                        TYPE *prev_x
                                        ) const
{
    if (number_of_rows != number_of_columns)
        throw std::exception("dense_matrix<TYPE>::GS_iteration: matrix is not squered");
    if (!x || !b || !prev_x)
        throw std::exception("dense_matrix<TYPE>::GS_iteration: invalid parameter");

    for (size_t row = 0; row < order; ++row)
    {
        x[row] = b[row];
        for (size_t col = 0; col < row; ++col)
            x[row] -= (a[row][col] * x[col]);
        for (size_t col = row + 1; col < order; ++col)
            x[row] -= (a[row][col] * prev_x[col]);
        x[row] *= omega;
        x[row] /= a[row][row];

        x[row] += (1 - omega) * prev_x[row];
    }
}
//------------------------------------------------------------------------------------------------ iterative_refinement
/**
*  Function improves solution of the equation Ax = b using specified method
*
*  @param method               - [in] describes a method by which we will approximate the exact solution
*  @param x                    - [in/out] searching solution (\in TYPE^number_of_columns)
*  @param b                    - [in] right sight vector of solving equation (\in TYPE^number_of_rows)
*  @param acc                  - [in] required accuracy (will not always be achieved)
*  @param max_it               - [in] maximum number of allowed iterations
*  @param DS  (default = NULL) - [in] optional, the input matrix A (before any changes like GE or something),
*                                it's used for calculating the residual vector in case main object of matrix
*                                is modified by GE, HOUSEHOLDER or something
*  @param d - (default = NULL) - [in] optional helpfull vector used to counting next iteration
*                                allocation and transmision of this parameter is recommended for large arrays,
*                                if y is not transmitted then function allocates it and deletes it itself
*  @param r - (default = NULL) - [in] optional helpfull residual vector used to counting next iteration
*                                allocation and transmision of this parameter is recommended for large arrays,
*                                if y is not transmitted then function allocates it and deletes it itself
*  @param y - (default = NULL) - [in] optional helpfull vector used to transmitted to function solve_LU
*                                allocation and transmision of this parameter is recommended for large arrays,
*                                if y is not transmitted then function allocates it and deletes it itself
*
*  @throw exception            - when matrix not square
*/
//---------------------------------------------------------------------------------------------------------------------
template <typename TYPE>
void dense_matrix<TYPE>::iterative_refinement( const EQUATION_METHOD method,
                                               TYPE *x,
                                               TYPE *b,
                                               double acc,
                                               size_t max_it,
                                               const dense_matrix<TYPE>* DS,
                                               TYPE *d,
                                               TYPE *r,
                                               TYPE *y
                                               )const
{
    if (number_of_columns != number_of_rows)
        throw std::exception("iterative_refinement fail: Matrix is not squared");

    const size_t N = number_of_columns;

    bool d_alloc = false; if (d == NULL) { d = new TYPE[N]; d_alloc = true; }
    bool r_alloc = false; if (r == NULL) { r = new TYPE[N]; r_alloc = true; }
    bool y_alloc = false; if (y == NULL) { y = new TYPE[N]; y_alloc = true; }

    if (DS != NULL)
        DS->count_rasidual_vector(x, b, r);
    else
        this->count_rasidual_vector(x, b, r);

    double v_norm = vector_norm(r, N);
    double new_v_norm;
    size_t iteration = 0;

    // int while condition are contained 2 conditins to stop the calculations,
    // third condition is implemented inside the loop
    // =======================================================================
    while (iteration < max_it && v_norm > acc)
    {
        switch (method)
        {
        case J:
            J_iteration(d, b, x);
            break;

        case GS:
            GS_iteration(d, b, x);
            break;

        case SOR:
            SOR_iteration(d, b, x);
            break;
        }

        if (DS != NULL)
            DS->count_rasidual_vector(d, b, r);
        else
            this->count_rasidual_vector(d, b, r);
        new_v_norm = vector_norm(r, N);

        // if norm of new residual vector is less then previous then accept new solution
        // =============================================================================
        if (new_v_norm < v_norm)
        {
            for (size_t i = 0; i < N; ++i)
                x[i] = d[i];
            v_norm = new_v_norm;
            iteration++;
        }
        // otherwise keep previous solution
        // ================================
        else
            break;
    }

    if (d_alloc) { delete[] d; d = NULL; }
    if (r_alloc) { delete[] r; r = NULL; }
    if (y_alloc) { delete[] y; y = NULL; }
}

//--------------------------------------------------------------------------- load_scheme_from_file
/**
*  Method used to loading matrix from file.
*
*  @param file_name         - [in] file name to read data
*/
/*  example format of input file:
*   _______________
*  |5              |  <-  number of rows of loading matrix
*  |5              |  <-  number of columns of loading matrix
*  |1  3  0  0  0  |
*  |0  0  2  0  1  |
*  |0  1  0  5  0  |
*  |0  0  0  1  0  |
*  |0  0  1  0  2  |
*  |_______________|
*
*  @throw exception     - if specified file not found
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE2>
void load_matrix_from_file(const char* file_name, dense_matrix<TYPE2>* DM)
{
}
template <>
void load_matrix_from_file(const char* file_name, dense_matrix<double>* DM)
{
char buffer[128];
std::ifstream input_file;
size_t n_of_rows, n_of_cols;


    input_file.open(file_name);

    if (input_file.is_open())
    {
        // First of all get sizes of matrix
        // ================================
        input_file.getline(buffer, 128);
        n_of_rows = atoi(buffer);
        input_file.getline(buffer,128);
        n_of_cols = atoi(buffer);

        DM->resize(n_of_rows, n_of_cols);

        // Read all elements in order
        // ==========================
        int row = 0, col = 0;
        while (!input_file.eof() && row != n_of_rows)
        {
            input_file >> buffer;
            DM->set_value(atof(buffer), row, col);
            ++col;
            if (col == n_of_cols)
            {
                col = 0;
                ++row;
            }
        }
        input_file.close();
    }
    else
    {
        std::string error("load_scheme_from_file fail: can not open the file: ");
        error += std::string(file_name);
        throw std::exception(error.c_str());
    }

}

//------------------------------------------------------------------------------------ dense_matrix
/// standard outstream operator
///
/// @param out            - out stream
/// @param DM             - input dense matrix
///
//-------------------------------------------------------------------------------------------------
template <typename TYPE2>
std::ostream& operator<< (std::ostream& out, const dense_matrix<TYPE2>& DM)
{
const size_t manip_double = 10;

    out << DM.number_of_rows << std::endl;
    out << DM.number_of_columns << std::endl;
    for (size_t row = 0; row < DM.number_of_rows; ++row)
    {
        for (size_t col = 0; col < DM.number_of_columns; ++col)
            out << std::setw(manip_double) << DM.a[row][col];
        out << std::endl;
    }
return out;
}

}; // std_math
#endif

