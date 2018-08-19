#ifndef MATRIX_STORAGE_SCHEMES
#define MATRIX_STORAGE_SCHEMES

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <algorithm>


//===============================================================================================//
/**
* \details this file contains very usefull classes to handle with sparse matrices
* \author £ukasz Jastrz¹b, email: luk.jastrzab@gmail.com
* \date 2012
* \version 1.0
* \copyright GNU Public Licence
*/
//===============================================================================================//


template <typename T>
double absolute_value(T number)
{
    return fabs(number);
}
template <>
double absolute_value(std::complex<float> number)
{
    return sqrt(number.real() * number.real() + number.imag() * number.imag());
}
template <>
double absolute_value(std::complex<double> number)
{
    return sqrt(number.real() * number.real() + number.imag() * number.imag());
}

//--------------------------------------------------------------------------- STANDARD VECTOR NORMS
// Euklides vector norm
// ====================
template <typename TYPE>
double vector_norm( TYPE *v,
                    size_t v_size
                   )
{
    double result = 0;
    for (size_t i = 0; i < v_size; ++i)
        result += (v[i] * v[i]);
    return sqrt(result);
}
template <>
double vector_norm( std::complex<float> *v,
                    size_t v_size
                    )
{
    double result = 0;
    for (size_t i = 0; i < v_size; ++i)
        result += v[i].real() * v[i].real() + v[i].imag() * v[i].imag();
    return sqrt(result);
}
template <>
double vector_norm( std::complex<double> *v,
                    size_t v_size
                    )
{
    double result = 0;
    for (size_t i = 0; i < v_size; ++i)
        result += v[i].real() * v[i].real() + v[i].imag() * v[i].imag();
    return sqrt(result);
}

/*************************************************************************************************/
/*                                                                                               */
/*                                     INPUT STORAGE SCHEME                                      */
/*                                                                                               */
/*************************************************************************************************/

template <typename TYPE>
class input_storage_scheme
{
public:
    /// number of added elements - sizes of arrays: AORIG, RNORIG and CNORIG
    const size_t NNZ;
    /// sizes of stored matrix
    const size_t number_of_rows, number_of_columns;
    /// not always real order of a matrix, just min(number_of_rows,number_of_columns)
    const size_t order;

private:
    /// array of ORIGinal values of the input matrix A
    std::vector<TYPE> AORIG;
    /// array of ORIGinal row numbers of the input matrix A (indexed from 0)
    std::vector<int> RNORIG;
    /// array of ORIGinal column numbers of the input matrix A (indexed from 0)
    std::vector<int> CNORIG;

public:
    /// default constructor
    input_storage_scheme()
        :
        NNZ(0),
        order(0),
        number_of_rows(0),
        number_of_columns(0)
    {
    }
    /// copy constructor
    input_storage_scheme(input_storage_scheme<TYPE>& ISS);
    /// constructor which requires sizes of matrix
    input_storage_scheme(size_t _number_of_rows, size_t _number_of_columns)
        :
        number_of_rows(_number_of_rows),
        number_of_columns(_number_of_columns),
        order(_number_of_rows < _number_of_columns ? _number_of_rows : _number_of_columns),
        NNZ(0)
    {
    }
    /// method to adding elements
    void add_element(TYPE value, size_t row, size_t col);

    /// method counts residual vector r_i = Ax_i - b
    void count_rasidual_vector(TYPE *x, TYPE *b, TYPE *r) const;
    /// assigned operator
    input_storage_scheme<TYPE>& operator= (input_storage_scheme<TYPE>& ISS);
private:

    /// Definition of basic out_stream operator
    template <typename TYPE2>
    friend std::ostream& operator<< (std::ostream& out, const input_storage_scheme<TYPE2>& ISS);
    /// Method uploading matrix from the file
    template <typename TYPE2>
    friend void load_matrix_from_file(const char* file_name, input_storage_scheme<TYPE2>* ISS);
    /// Method uploading scheme from the file
    template <typename TYPE2>
    friend void load_scheme_from_file(const char* file_name, input_storage_scheme<TYPE2>* ISS);
    /// Friend class basing on input_storage_scheme
    template <typename TYPE2>
    friend class dynamic_storage_scheme;

    /// External test functions
    template <typename TYPE2>
    friend void permute_input_matrix_elements_test(input_storage_scheme<TYPE> *ISS);
    /// 
    template <typename TYPE2>
    friend void CheckAxb(const input_storage_scheme<TYPE2> &ISS, const TYPE2 *x, TYPE2 *b);

};
//---------------------------------------------------------------------------- input_storage_scheme
/**
*  Copy constructor
*
*  @param ISS            - [in] object to be copied
*
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE> 
input_storage_scheme<TYPE>::input_storage_scheme(input_storage_scheme<TYPE>& ISS)
:
number_of_rows(ISS.number_of_rows),
number_of_columns(ISS.number_of_columns),
order(ISS.order),
NNZ(ISS.NNZ)
{
    AORIG = ISS.AORIG;
    RNORIG = ISS.RNORIG;
    CNORIG = ISS.CNORIG;
}
//------------------------------------------------------------------------------------- add_element
/**
*  Method used to adding elements
*
*  @param value            - [in] main content of added element
*  @param row              - [in] row number of added element
*  @param col              - [in] column number of added element
*
*  @throw out_of_range
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE> 
void input_storage_scheme<TYPE>::add_element( TYPE value,
                                              size_t row,
                                              size_t col
                                              )
{
    if (row >= number_of_rows || col >= number_of_columns)
        throw std::out_of_range("input_storage_scheme<TYPE>::add_element: wrong row/col number");

    AORIG.push_back(value);
    RNORIG.push_back(row);
    CNORIG.push_back(col);
    const_cast<size_t&>(NNZ)++;
}

//--------------------------------------------------------------------------- count_rasidual_vector
/**
*  Method used to counting residual vector
*
*  @param x                - [in]  obtained solution
*  @param b                - [in]  right sight vestor of equation Ax=b
*  @param r                - [out] counting residual vestor
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void input_storage_scheme<TYPE>::count_rasidual_vector( TYPE *x,
                                                        TYPE *b,
                                                        TYPE *r
                                                        ) const
{
    for (size_t row = 0; row < number_of_rows; ++row)
        r[row] = -b[row];
    for (size_t idx = 0; idx < NNZ; ++idx)
        r[RNORIG[idx]] += (x[CNORIG[idx]] * AORIG[idx]);
}

//--------------------------------------------------------------------------- count_rasidual_vector
/**
*  The assignment operator
*
*  @param ISS                - [in] the assigned input scheme
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
input_storage_scheme<TYPE>& input_storage_scheme<TYPE>:: operator= (input_storage_scheme<TYPE>& ISS)
{
    const_cast<size_t&>(number_of_rows) = ISS.number_of_rows;
    const_cast<size_t&>(number_of_columns) = ISS.number_of_columns;
    const_cast<size_t&>(order) = ISS.order;
    const_cast<size_t&>(NNZ) = ISS.NNZ;
    AORIG = ISS.AORIG;
    RNORIG = ISS.RNORIG;
    CNORIG = ISS.CNORIG;
    return *this;
}
//-------------------------------------------------------------------------------------- operator<<
/**
*  Standard outstream operator
*
*  @param out                           -[in/out] standard outstream operator
*  @param ISS                           -[in] input storage scheme
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
std::ostream& operator<< ( std::ostream& out,
                           const input_storage_scheme<TYPE>& ISS
                           )
{
const int manip_typ = 6;
const int manip_int = 3;

    out << ISS.number_of_rows << std::endl
        << ISS.number_of_columns << std::endl;

    for (size_t i = 0; i < ISS.AORIG.size(); i++)
        out << std::setw(manip_typ) << ISS.AORIG[i] << std::setw(manip_int) << ISS.RNORIG[i]
            << std::setw(manip_int) << ISS.CNORIG[i] << std::endl;

    return out;
}
//--------------------------------------------------------------------------- load_matrix_from_file
/**
*  Method used to loading matrix from file.
*  In file must be saved also zero elements, functions avoid them during adding elements
*
*  @param file_name         - [in]  file name to read data
*  @param ISS               - [out] memorized input storage scheme
*
*  example format of input file:
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
/// Declaration
template <typename TYPE>
void load_matrix_from_file( const char* file_name,
                            input_storage_scheme<TYPE>* ISS
                            );
/// template for specjalization in loading elements of specified type
template <typename TYPE>
void load_matrix_from_file( const char *file_name,
                            input_storage_scheme<TYPE>* ISS
                            )
{
std::ifstream input_file;
char buffer[128];

    // Open input file
    // ===============
    input_file.open(file_name);

    if (input_file.is_open())
    {
        // First of all get sizes of matrix
        // ================================
        input_file.getline(buffer, 128);
        const_cast<size_t&>(ISS.number_of_rows) = atoi(buffer);
        input_file.getline(buffer,128);
        const_cast<size_t&>(ISS.number_of_columns) = atoi(buffer);
        ISS.order = (ISS.number_of_rows < ISS.number_of_columns ? ISS.number_of_rows : ISS.number_of_columns);

        // ======================================================= //
        // Second - ADD YOUR OWN CODE TO STORING SPECIFIED OBJECTS //
        // ======================================================= //
        input_file.close();

        // =======================================================
        throw std::exception("load_matrix_from_file fail: not allowed specialized function");
        // becouse the specialization should be implemented in this particular moment :]
    }
    else
    {
        std::string error("load_matrix_from_file fail: can not open the file: ");
        error += std::string(file_name);
        throw std::exception(error.c_str());
    }

}
///  Specjalization for double
template <>
void load_matrix_from_file( const char *file_name,
                            input_storage_scheme<double>* ISS
                            )
{
std::ifstream input_file;
char buffer[128];

    // Open input file
    // ===============
    input_file.open(file_name);

    if (!input_file.is_open())
    {
        input_file.close();
        std::string error("load_matrix_from_file fail: can not open the file ");
        error += std::string(file_name);
        throw std::exception(error.c_str());
    }

    // First of all get sizes of matrix
    // ================================
    input_file.getline(buffer, 128);
    const_cast<size_t&>(ISS->number_of_rows) = atoi(buffer);
    input_file.getline(buffer,128);
    const_cast<size_t&>(ISS->number_of_columns) = atoi(buffer);
    const_cast<size_t&>(ISS->order) = (ISS->number_of_rows < ISS->number_of_columns ? ISS->number_of_rows : ISS->number_of_columns);

    // Read all elements in order
    // ==========================
    int row = 0, col = 0;
    while (!input_file.eof() && row != ISS->number_of_rows)
    {
        input_file >> buffer;
        const double val = atof(buffer);
        if (val != 0)
           ISS->add_element(val, row, col);

        col++;
        if (col == ISS->number_of_columns)
        {
            col = 0;
            row ++;
        }
    }

input_file.close();
}
//--------------------------------------------------------------------------- load_matrix_from_file
/**
*  Method used to loading scheme from file.
*  In file must be saved also non zero elements in a simple list
*
*  @param file_name         - [in]  file name to read data
*  @param ISS               - [out] memorized input storage scheme
*
*  example format of input file:
*   _______________
*  |5              |  <-  number of rows of loading matrix
*  |5              |  <-  number of columns of loading matrix
*  |1  0  0        |  <-  TYPE    row     col
*  |5  1  2        |  <-  TYPE    row     col
*  |3  1  0        |  <-  TYPE    row     col
*  |2  2  0        |  <-  TYPE    row     col
*  |5  3  1        |  <-  TYPE    row     col
*  |_______________|
*
*  @throw exception     - if specified file not found
*/
//-------------------------------------------------------------------------------------------------
/// Declaration
template <typename TYPE>
void load_scheme_from_file( const char* file_name,
                            input_storage_scheme<TYPE>* ISS
                            );
/// Declaration
template <typename TYPE>
void load_scheme_from_file( const char* file_name,
                            input_storage_scheme<TYPE>* ISS
                            )
{
std::ifstream input_file;
char buffer[128];

    // Open input file
    // ===============
    input_file.open(file_name);

    if (!input_file.is_open())
    {
        input_file.close();
        std::string error("load_scheme_from_file fail: can not open the file ");
        error += std::string(file_name);
        throw std::exception(error.c_str());
    }

    // First of all get sizes of matrix
    // ================================
    input_file.getline(buffer, 128);
    const_cast<size_t&>(ISS->number_of_rows) = atoi(buffer);
    input_file.getline(buffer,128);
    const_cast<size_t&>(ISS->number_of_columns) = atoi(buffer);
    const_cast<size_t&>(ISS->order) = (ISS->number_of_rows < ISS->number_of_columns ? ISS->number_of_rows : ISS->number_of_columns);

    // Read all elements in order
    // ==========================
    while (!input_file.eof())
    {
        input_file >> buffer;
        const double val = atof(buffer);
        input_file >> buffer;
        const int row = atoi(buffer);
        input_file >> buffer;
        const int col = atoi(buffer);
        if (val != 0)
           ISS->add_element(val, row, col);
    }

input_file.close();
}
//-------------------------------------------------------------------------------------------------
/**
*   Method random permuts elements of the input matrix
*
*  @param ISS -      input storage scheme to be tested
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void permute_input_matrix_elements_test(input_storage_scheme<TYPE> *ISS)
{
    srand((size_t)time(NULL));
    const size_t N = ISS->NNZ;
    TYPE2 val;
    int row, col;

    for (size_t idx1 = 0; idx1 < N; ++idx1)
        for (size_t idx2 = idx1 + 1; idx2 < N; ++idx2)
        {
            if (rand() % 2 == 0)
            {
                val = ISS->AORIG[idx1];
                row = ISS->RNORIG[idx1];
                col = ISS->CNORIG[idx1];
                ISS->AORIG[idx1] = ISS->AORIG[idx2];
                ISS->RNORIG[idx1] = ISS->RNORIG[idx2];
                ISS->CNORIG[idx1] = ISS->CNORIG[idx2];
                ISS->AORIG[idx2] = val;
                ISS->RNORIG[idx2] = row;
                ISS->CNORIG[idx2] = col;
            }
        }
}
//-------------------------------------------------------------------------------------------------
/**
*   Method checks soltion x - (it multiplies A * x and returns result b)
*
*  @param ISS -      [in] input storage scheme to be tested
*  @param x -        [in] solution x
*  @param b -        [out] result
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE2>
void CheckAxb(const input_storage_scheme<TYPE2> &ISS, const TYPE2 *x, TYPE2 *b)
{
       for (size_t i = 0; i < ISS.number_of_rows; ++i)
           b[i] = 0;
       for (size_t i = 0; i < ISS.NNZ; ++i)
           b[ISS.RNORIG[i]] += ISS.AORIG[i] * x[ISS.CNORIG[i]];
}


//------------------------- ENUMERATES USING ONLY BY DYNAMIC STORAGE SCHEME -----------------------
// Type definition for state of dynamic_storage_scheme
// ===================================================
typedef
enum _DYNAMIC_STATE
{
    ROL_INIT,
    COL_INIT,
    ITERATIVE,
    LU_DECOMPOSED,
    QR_DECOMPOSED
}
DYNAMIC_STATE;

// Type definition for using pivotal strategy
// ==========================================
typedef
enum _PIVOTAL_STRATEGY
{
    ONE_ROW_SEARCHING,
    MARKOWITZ_COST,
    FILLIN_MINIMALIZATION
}
PIVOTAL_STRATEGY;

// Enumerates used for check if storing fillins was successful
// ===========================================================
enum STORING_STATUS
{
    STORING_SUCCESS,
    STORING_FAIL
};

// Enums for marking free/not free positions and empty row/column packages
// =======================================================================
enum MARKER
{
    EMPTY = -2,
    FREE = -1,
    NFREE = 1
};

/*************************************************************************************************/
/*                                                                                               */
/*                                   DYNAMIC STORAGE SCHEME                                      */
/*                                                                                               */
/*************************************************************************************************/
template <typename TYPE>
class dynamic_storage_scheme
{
private:
    //======== MATRIX - BASIC INFORMATIONS ========
    /// sizes of stored matrix
    const size_t number_of_rows, number_of_columns;
    /// mostly  = min(number_of_rows, number_of_columns)
    const size_t order;

    //============== Row-Ordered List (ROL) =============
    size_t NROL;    /// size of row-ordered list (not number of stored elements)
    size_t LROL;    /// last not free position in ROL
    size_t CROL;    /// number of non-zeros actualy stored in ROL (calculated control)
    TYPE *ALU;      /// array of values of the elements stored in scheme (indexed from 0)
    int *CNLU;      /// array of column numbers of the elements stored in scheme (indexed from 0)

    //=============== Column-Ordered List ===============
    size_t NCOL;    /// size of column-ordered list (not number of stored elements)
    size_t LCOL;    /// last not free position in COL
    size_t CCOL;    /// number of non-zeros actualy stored in COL (calculated control)
    int *RNLU;      /// array of row numbers of the elements stored in scheme (indexed from 0)

    //============== INTEGRITY ARRAYS =============
    size_t NHA;     /// size of integrity arrays
    TYPE *PIVOT;    /// table of pivots

    int **HA; /// Array of pointers and permutations
    /**
        HA[i][0] - place to store working indexes during algorythmization

        HA[r][1] - first element in r-th row in ROL
        HA[r][2] - midle pointer in r-th row in ROL
        HA[r][3] - last element in r-th row in ROL
        HA[c][4] - first element in c-th column in COL
        HA[c][5] - midle pointer in c-th column in COL
        HA[c][6] - last element in c-th column in COL

        // permutations part
        HA[r][7] - original number of row, which is placed on r-th position now
        HA[r][8] - number of row, where is placed r-th original row
        HA[c][9] - original number of column, which is placed on c-th position now
        HA[c][10] - number of column, where is placed c-th original column
    */

    DYNAMIC_STATE dynamic_state;    /// variable indicating current state of the scheme
    double omega;                   /// relaxation parameter for SOR method
    std::string logged_errors;      /// output string log all errors and warnings obtained during the computation


public:

    /// Constructor - input_storage_scheme require and two floats that determine sizes of storage lists
    ///               in flollowing way: NROL = mult1 * ISS->NNZ, NCOL = mult2 * NROL
    dynamic_storage_scheme(const input_storage_scheme<TYPE>& ISS, double mult1, double mult2 = 0.7, DYNAMIC_STATE _dynamic_state = ROL_INIT);
    /// Destructor
    ~dynamic_storage_scheme();

    /// Method gets string with logged warnings and errors
    std::string get_logged_errors(void) { return logged_errors; }
    /// Method clears logged_errors
    void clear_logged_errors(void) { logged_errors.clear(); }
    /// Method to dump the contents of the schema
    void print_scheme_to_file(const char *file_name);
    /// Method used for decomposition of matrix (Gauss elimination)
    void LU_decomposition(PIVOTAL_STRATEGY strategy, size_t _search, double _mult, double eps, bool pre_sort = false);
    /// Method solves LU problem (LU_decomposition is needed to call before)
    void solve_LU(TYPE *x, TYPE *b, TYPE *y = NULL)const;
    /// Method improves the accuracy of the solution
    void iterative_refinement(const input_storage_scheme<TYPE>& ISS, TYPE *x, TYPE *b, double acc, size_t max_it) const;
    /// Method prepares matrix to SOR iterations
    void iterative_preparation(void);
    /// Method sets relaxation parameter omega
    void set_omega(double _omega)
    {
        omega = _omega;
    }
    /// Method performs one iteration of SOR method
    void SOR_iteration(TYPE *x, TYPE *b, TYPE *prev_x) const;


    /// Method returns ostream with sparsity patern
    void print_sparsity_pattern(const char *file_name);

private:
    // Disable default constructor
    dynamic_storage_scheme();

    //============== DYNAMIC SCHEME MANIPULATORS ==============
    /// Function permuts row lying on pos1 position with row lying on pos2 position
    void permute_rows(size_t pos1, size_t pos2);
    /// Function permuts column lying on pos1 position with column lying on pos2 position
    void permute_columns(size_t pos1, size_t pos2);
    /// Function use to storing fillins in row ordered list, params row and col are current position in matrix
    int store_fillin_ROL(TYPE val, int row, int col, bool garbage_on = true);
    /// Function performs the organize of the elements in a compact structure in ROL
    void garbage_collection_in_ROL(void);
    /// Function use to storing fillins in column ordered list, params row and col are current position in matrix
    int store_fillin_COL(int row, int col, bool garbage_on = true);
    /// Function performs the organize of the elements in a compact structure in COL
    void garbage_collection_in_COL(void);
    /// Method brings the choosen element on choosen row to begin of active part of it and incress active begine pointer (so element is inactive)
    void deactivate_element_in_ROL(size_t index, size_t row);
    /// Method switch last element with the specified and sets FREE to cuurent last and decrease end row pointer
    void eliminate_element_in_ROL(size_t index, size_t row);
    /// Method switch last element with the specified and sets FREE to cuurent last and decrease end row pointer
    void eliminate_element_in_COL(size_t index, size_t col);
    /// Method deletes choosen column package in column ordered list
    void delete_column_package_in_COL(size_t col);


    //============== PIVOTAL STRATEGIES ==============
    /// Strategy choosing pivot in s-th row only (pre-rows sort is aplied before)
    void choose_pivot_by_ONE_ROW_SEARCHING(size_t stage, double mult);
    /// Strategy based on Markowitz cost function
    void choose_pivot_by_MARKOWITZ_COST(size_t stage, size_t search, double mult);
    /// Strategy based on local minimalization of arising fillins
    void choose_pivot_by_FILLIN_MINIMALIZATION(size_t stage, size_t search, double mult);

    /// Quick sort method for sorting rows
    void sort_rows(int l = 0, int r = 0);
    /// Method counts the fillin cost in assumption that specified element was choosen as a pivot
    int count_fillin_cost(size_t index, size_t row);

    //=============== FRIEND FUNCTIONS ===============
    /// Definition of basic out_stream operator
    template <typename TYPE2>
    friend std::ostream& operator<<(std::ostream& out, const dynamic_storage_scheme<TYPE2>& DSS);


    /// ==================== FUNCTIONS WRITTEN ONLY FOR TESTS ====================
    /// Method checks if there is not any mistakes in scheme structure
    int check_integrity_test(void)const;

};


//-------------------------------------------------------------------------- dynamic_storage_scheme
/**
*  The only one constructor which requires following parameters:
*
*  @param ISS              - [in] The input scheme under which creates a dynamic scheme
*  @param mult1            - [in] sizeof(ROL) = mult1 * sizeof(ISS)
*  @param mult2 = 0.7      - [in] sizeof(COL) = mult2 * sizeof(ROL)
*
*  @throw bad_alloc        - when scheme couldn't be allocate
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
dynamic_storage_scheme<TYPE>::
dynamic_storage_scheme( const input_storage_scheme<TYPE>& ISS,
                                                      double mult1,
                                                      double mult2,
                                                      DYNAMIC_STATE _dynamic_state
                                                      )
:
number_of_rows(ISS.number_of_rows),
number_of_columns(ISS.number_of_columns),
order(ISS.number_of_rows < ISS.number_of_columns ? ISS.number_of_rows : ISS.number_of_columns),
ALU(NULL),
CNLU(NULL),
RNLU(NULL),
PIVOT(NULL),
HA(NULL),
omega(1),
dynamic_state(_dynamic_state)
{
    if (dynamic_state != ROL_INIT && dynamic_state != COL_INIT)
        dynamic_state = ROL_INIT;

    const size_t  NNZ = ISS.NNZ;

    // Allocate additional tabs for check the content each row/col - pack in COL/ROL
    // =============================================================================
    bool *rows_content_check = NULL,
         *cols_content_check = NULL;

    // set sizes of allocated memory for state ROL_INIT
    // ================================================
    if (dynamic_state == ROL_INIT)
    {
        NROL = static_cast<size_t>(mult1 * NNZ);
        if (NROL < NNZ) NROL = NNZ;
        NCOL = static_cast<size_t>(mult2 * NROL);
        if (NCOL < NNZ) NCOL = NNZ;
    }
    // set sizes of allocated memory for state COL_INIT
    // ================================================
    else
    {
        NCOL = static_cast<size_t>(mult1 * NNZ);
        if (NCOL < NNZ) NCOL = NNZ;
        NROL = static_cast<size_t>(mult2 * NCOL);
        if (NROL < NNZ) NROL = NNZ;
    }
    NHA = (number_of_rows > number_of_columns ? number_of_rows : number_of_columns);

    // Allocation of memory
    // ====================
    try
    {
        if (dynamic_state == ROL_INIT)
            ALU = new TYPE[NROL];
        else
            ALU = new TYPE[NCOL];
        CNLU = new int[NROL];
        RNLU = new int[NCOL];
        PIVOT = new TYPE[NHA];
        HA = new int*[NHA];
        for (size_t i = 0; i < NHA; ++i)
            HA[i] = NULL;
        for (size_t i = 0; i < NHA; ++i)
            HA[i] = new int[11];

        // temporary additional fields
        // ===========================
        rows_content_check = new bool[number_of_rows];
        cols_content_check = new bool[number_of_columns];
    }
    // Release everything if exception ocured
    // ======================================
    catch (std::bad_alloc& e)
    {
        delete[] ALU;
        delete[] CNLU;
        delete[] RNLU;
        delete[] PIVOT;
        if (HA != NULL)
            for (size_t i = 0; i < NHA; i++)
                delete[] HA[i];
        delete[] HA;
        delete[] rows_content_check;
        delete[] cols_content_check;
        throw e;
    }
    // ====================================
    // Few steps of optimal initialization 
    // ====================================

    // Second step - count number of elements in each row/column package
    //==================================================================
    for (size_t i = 0; i < NHA; ++i)
        HA[i][3] = HA[i][6] = 0;
    for (size_t i = 0; i < NNZ; ++i)
    {
        HA[ISS.RNORIG[i]][3]++;
        HA[ISS.CNORIG[i]][6]++;
    }

    // Set temporary initializations for poroper allocation of memory
    // in case there are some empty row/column packages (matrix will be singular!)
    // ===========================================================================
    for (size_t row = 0; row < number_of_rows; row++)
        rows_content_check[row] = (HA[row][3] > 0 ? true : false);
    for (size_t col = 0; col < number_of_columns; col++)
        cols_content_check[col] = (HA[col][6] > 0 ? true : false);

    // Third step - setting row/column packages beginnings
    //====================================================
    HA[0][1] = HA[0][4] = 0;
    for (size_t i = 0; i < NHA - 1; ++i)
    {
        HA[i+1][1] = HA[i][1] + HA[i][3];
        HA[i+1][4] = HA[i][4] + HA[i][6];
        HA[i][3] = HA[i][6] = 0;
    }
    HA[NHA - 1][3] = HA[NHA - 1][6] = 0;

    // initialization of ROL and COL in ROL_INIT dynamic_state
    // =======================================================
    if (dynamic_state == ROL_INIT)
    {
        // Fourth step - storing elements in ROL
        //======================================
        for (size_t i = 0; i < NNZ; ++i)
        {
            const int k = ISS.RNORIG[i],
                      l = HA[k][1] + HA[k][3];
            CNLU[l] = ISS.CNORIG[i];
            ALU[l] = ISS.AORIG[i];
            ++HA[k][3];
        }
        // Fifth step - setting middle and end pointers for ROL and initiation of COL
        //===========================================================================
        for (size_t i = 0; i < NHA; ++i)
        {
            HA[i][2] = HA[i][1];
            HA[i][3] = HA[i][1] + HA[i][3] - 1;
            for (int j = HA[i][1]; j <= HA[i][3]; j++)
            {
                const int l = CNLU[j];
                RNLU[HA[l][4] + HA[l][6]] = i;
                ++HA[l][6];
            }
        }
        // Sixth step - setting middle and end pointers for COL
        //=====================================================
        for (size_t i = 0; i < NHA; ++i)
        {
            HA[i][5] = HA[i][4];
            HA[i][6] = HA[i][4] + HA[i][6] - 1;
        }
    }
    // initialization of ROL and COL in COL_INIT dynamic_state
    // =======================================================
    else
    {
        // Fourth step - storing elements in COL
        //======================================
        for (size_t i = 0; i < NNZ; ++i)
        {
            const int k = ISS.CNORIG[i],
                      l = HA[k][4] + HA[k][6];
            RNLU[l] = ISS.RNORIG[i];
            ALU[l] = ISS.AORIG[i];
            ++HA[k][6];
        }
        // Fifth step - setting middle and end pointers for ROL and initiation of COL
        //===========================================================================
        for (size_t i = 0; i < NHA; ++i)
        {
            HA[i][5] = HA[i][4];
            HA[i][6] = HA[i][4] + HA[i][6] - 1;
            for (int j = HA[i][4]; j <= HA[i][6]; ++j)
            {
                const int l = RNLU[j];
                CNLU[HA[l][1] + HA[l][3]] = i;
                ++HA[l][3];
            }
        }
        // Sixth step - setting middle and end pointers for COL
        //=====================================================
        for (size_t i = 0; i < NHA; ++i)
        {
            HA[i][2] = HA[i][1];
            HA[i][3] = HA[i][1] + HA[i][3] - 1;
        }
    }

    // Seventh step - mark free positions in allocated lists
    // settin ALU[i] = 0 for i-th free position is unnecessary
    //========================================================
    for (size_t i = NNZ; i < NROL; ++i)
        CNLU[i] = FREE;
    for (size_t i = NNZ; i < NCOL; ++i)
        RNLU[i] = FREE;

    // Seting permutations to identical
    // ================================
    for (size_t i = 0; i < NHA; ++i)
        for (size_t j = 7; j < 11; ++j)
            HA[i][j] = i;

    // clear operating memory
    // ======================
    for (size_t i = 0; i < NHA; ++i)
    {
        HA[i][0] = FREE;
        PIVOT[i] = 0;
    }

    // Setting last not-free markers and counters for both lists
    // =========================================================
    LROL = LCOL = NNZ - 1;
    CROL = CCOL = NNZ;

    // Use of content checks for set empty rows/cols
    // =============================================
    for (size_t row = 0; row < number_of_rows; ++row)
        if (!rows_content_check[row])
        {
            HA[row][1] = HA[row][2] = FREE;
            HA[row][3] =  EMPTY;
        }
    for (size_t col = 0; col < number_of_columns; ++col)
        if (!cols_content_check[col])
        {
            HA[col][4] = HA[col][5] = FREE;
            HA[col][6] = EMPTY;
        }

    // delete redundant data
    // =====================
    delete[] rows_content_check;
    delete[] cols_content_check;
}

//------------------------------------------------------------------------- ~dynamic_storage_scheme
/**
*  Deleting dynamic allocated (by operator new[]) memory
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
dynamic_storage_scheme<TYPE>::~dynamic_storage_scheme(void)
{
    delete[] ALU;
    delete[] CNLU;
    delete[] RNLU;
    delete[] PIVOT;
    if (HA != NULL)
        for (size_t i = 0; i < NHA; i++)
            delete[] HA[i];
    delete[] HA;
}

//---------------------------------------------------------------------------- print_scheme_to_file
/**
*  Method prints scheme to file
*  @param file_name         - [in] name of file to which scheme should be printed
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::print_scheme_to_file(const char *file_name)
{
std::ofstream outFile;

    // Try to open/create the file
    // ===========================
    try
    {
        outFile.open(file_name);
    }
    catch(std::exception)
    {
        logged_errors += "dynamic_storage_scheme<TYPE>::print_scheme_to_file: can not open the file: ";
        logged_errors += std::string(file_name);
        logged_errors += "\n";
        outFile.close();
        return;
    }
    outFile << *this;

outFile.close();
}



//-------------------------------------------------------------------------------- LU_decomposition
/**
*  Main functionality of the dynamic_storage_scheme, decompusition perfomed by Gauss
*  elimination algorithm.
*
*  @param strategy                     - [in] parameter discribing which strategy should be
*                                        used to choise pivot
*  @param search                       - [in] number of rows to be searched during pivot choosing
*                                        (not valid for ONE_ROW_SEARCHING strategy)
*  @param mult                         - [in] stability parameter for choosing pivot
*  @param eps                          - [in] parameter used for deleting to small elements,
*                                        or not storing to small fillins, that meets (fabs < eps)
*  @param pre_sort - (default = false) - [in] flag indicating if pre row sorting should be performed,
*                                        default it will be performed only if ONE_ROW_SEARCHING
*                                        strategy is in use
*
*  @throw exception                    - when matrix is not square, singular or there is not
*                                        enough memory in ROL/COL
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::
LU_decomposition( PIVOTAL_STRATEGY strategy,
                  size_t _search,
                  double _mult,
                  double eps,
                  bool pre_sort
                 )
{
    if (dynamic_state != ROL_INIT)
        throw std::exception("dynamic_storage_scheme<TYPE>::LU_decomposition: ROL_INIT state is required\n");

    if (number_of_columns != number_of_rows)
        throw std::exception("dynamic_storage_scheme<TYPE>::LU_decomposition: matrix is not squared");

    const double mult = (_mult > 1 ? _mult : 1);
    const size_t search = (strategy == ONE_ROW_SEARCHING || _search < 1 ? 1 : _search);
    const int N = (int)number_of_rows - 1;

    if (pre_sort || strategy == ONE_ROW_SEARCHING)
        sort_rows();

    // main loop, over all stages of elimination
    // =========================================
    for (int stage = 0; stage < N; ++stage)
    {
        // choose pivot and permute rows and columns
        // =========================================
        switch (strategy)
        {
            case ONE_ROW_SEARCHING:
                choose_pivot_by_ONE_ROW_SEARCHING(stage, mult);
                break;
            case MARKOWITZ_COST:
                choose_pivot_by_MARKOWITZ_COST(stage, search, mult);
                break;
            case FILLIN_MINIMALIZATION:
                choose_pivot_by_FILLIN_MINIMALIZATION(stage, search, mult);
                break;
        }
        if (absolute_value(PIVOT[stage]) < eps)
            throw std::exception("dynamic_storage_scheme<TYPE>::LU_decomposition: obtained singular matrix");

        const int eliminating_row = HA[stage][7];

        // over active part of column which is placed on stage-th position
        // ===============================================================
        const size_t stage_col = HA[stage][9];
        int col_end = HA[stage_col][6];
        for (int col_idx = HA[stage_col][5]; col_idx <= col_end; ++col_idx)
        {
            if (HA[RNLU[col_idx]][8] <= stage)
                continue;
            else
            {
                const int eliminated_row = RNLU[col_idx];

                // store elements of eliminating row to PIVOT on its column numbers positions
                // ==========================================================================
                for (int idx = HA[eliminating_row][2]; idx <= HA[eliminating_row][3]; ++idx)
                    PIVOT[HA[CNLU[idx]][10]] = ALU[idx]; 

                // first find element in eliminated row which has the same column number as pivot and count eliminator
                // ===================================================================================================
                TYPE eliminator;
                for (int idx = HA[eliminated_row][2]; idx <= HA[eliminated_row][3]; ++idx)
                    if (CNLU[idx] == HA[stage][9])
                    {
                        eliminator = ALU[idx] = ALU[idx] / PIVOT[stage];
                        deactivate_element_in_ROL(idx, eliminated_row);
                        break;
                    }

                // second, modify existing in eliminated row elements with the same column number as non-zeros in eliminating row
                // ==============================================================================================================
                for (int idx = HA[eliminated_row][2]; idx <= HA[eliminated_row][3]; ++idx)
                {
                    const size_t IDX = HA[CNLU[idx]][10];
                    if (absolute_value(PIVOT[IDX]) != 0)
                    {
                        ALU[idx] = ALU[idx] - eliminator * PIVOT[IDX];
                        PIVOT[IDX] = 0;

                        // if modified element is to small then eliminate it from COL and ROL
                        // ==================================================================
                        if (absolute_value(ALU[idx]) <= eps)
                        {
                            const int eliminated_col = CNLU[idx];
                            for (int idxc = HA[eliminated_col][5]; idxc <= HA[eliminated_col][6]; ++idxc)
                                if (RNLU[idxc] == eliminated_row)
                                {
                                    eliminate_element_in_COL(idxc, eliminated_col);
                                    break;
                                }

                            eliminate_element_in_ROL(idx, eliminated_row);
                            idx--;
                        }
                    }
                }
                // last part - add fillins
                // =======================
                for (int idx = HA[eliminating_row][2]; idx <= HA[eliminating_row][3]; ++idx)
                {
                    const size_t col_number = CNLU[idx];
                    TYPE fillin = - eliminator * PIVOT[HA[col_number][10]];
                    if (absolute_value(fillin) > eps)
                    {
                        if (store_fillin_ROL(fillin, eliminated_row, col_number) == STORING_FAIL)
                            throw std::exception("dynamic_storage_scheme<TYPE>::LU_decomposition: not enough memory in ROL");
                        if (store_fillin_COL(eliminated_row, col_number) == STORING_FAIL)
                            throw std::exception("dynamic_storage_scheme<TYPE>::LU_decomposition: not enough memory in COL");
                        PIVOT[HA[col_number][10]] = 0;
                    }
                }
            }
            // if garbage collection in COL is performed, then reindex col_idx and change col_end
            // ==================================================================================
            if (HA[stage_col][6] != col_end)
            {
                col_idx -= (col_end - HA[stage_col][6]);
                col_end = HA[stage_col][6];
            }
        }

        // delete column on stage position (it is useless now)
        // ===================================================
        delete_column_package_in_COL(stage_col);

        // delete elements from COL which rows has been eliminators,
        // this is necessery for correct work of pivotal strategies
        // =========================================================
        const int row_end = HA[HA[stage][7]][3];
        for (int rol_idx = HA[HA[stage][7]][2]; rol_idx <= row_end; ++rol_idx)
        {
            const int col = CNLU[rol_idx];
            for (int col_idx = HA[col][5]; col_idx <= HA[col][6]; ++col_idx)
                if (HA[RNLU[col_idx]][8] <= stage)
                {
                    eliminate_element_in_COL(col_idx, col);
                    break;
                }
        }
    }

    // choose last pivot using the simplest strategy
    // =============================================
    choose_pivot_by_ONE_ROW_SEARCHING(N, mult);

    // if last pivot is 0, it meens that matrix became singular
    // ========================================================
    if (absolute_value(PIVOT[N]) == 0)
       throw std::exception("dynamic_storage_scheme<TYPE>::LU_decomposition: obtained singular matrix");

    dynamic_state = LU_DECOMPOSED;
}

//--------------------------------------------------------------------------------------------- solve_LU
/**
*  Function solves equation LUx = b, where L and U are respectively down and upper
*  triangular matrices obtained by LU_decomposition method from input matrix,
*  so there is a need to call this method before
*
*  @param x                     - [out] searching solution (\in TYPE^number_of_columns)
*  @param b                     - [in] right sight vector of solving equation (\in TYPE^number_of_rows)
*  @param y  - (default = NULL) - [in] optional helpfull vector used to counting solution in two stages
*                                 Ly = b and Ux = y, allocation and transmision of this parameter 
*                                 is recommended for large arrays, if y is not transmitted then
*                                 function allocates it and deletes by self
*
*  @throw exception             - when matrix is not decomponed or matrix is not square
*/
//------------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::solve_LU( TYPE *x,
                                             TYPE *b,
                                             TYPE *y
                                             ) const
{
    if (dynamic_state != LU_DECOMPOSED)
        throw std::exception("dynamic_storage_scheme<TYPE>::solve_LU: LU_decomposition is needed before");
    if (number_of_columns != number_of_rows)
        throw std::exception("dynamic_storage_scheme<TYPE>::solve_LU: matrix is not squared");

    const size_t N = number_of_columns;

    bool y_alloc = false; if (y == NULL) { y = new TYPE[N]; y_alloc = true; }

    // first solve the equation Ly = b
    // ===============================
    y[0] = b[HA[0][7]];
    for (size_t row = 1; row < N; row++)
    {
        const int origRow = HA[row][7];
        y[row] = b[origRow];
        TYPE sum = 0;
        for (int idx = HA[origRow][1]; idx < HA[origRow][2]; idx++)
            sum += ALU[idx] * y[HA[CNLU[idx]][10]];
        y[row] -= sum;
    }
    // second solve the equation Ux = y
    // ================================
    x[HA[N - 1][9]] = y[N - 1] / PIVOT[N - 1];
    for (int row = N - 2; row >= 0; row--)
    {
        const int col = HA[row][9];
        const int origRow = HA[row][7];
        x[col] = y[row];
        for (int idx = HA[origRow][2]; idx <= HA[origRow][3]; idx++)
            x[col] -= ALU[idx] * x[CNLU[idx]];
        x[col] /= PIVOT[row];
    }

    if (y_alloc) delete[] y;
}

//------------------------------------------------------------------------------------------------ iterative_refinement
/**
*  Function improves solution of the equation Ax = b using specified method
*
*  @param ISS                  - [in] the input scheme under which the dynamic scheme has been created
*                                for the benefit of which the method is called
*  @param x                    - [in/out] searching solution (\in TYPE^number_of_columns)
*  @param b                    - [in] right sight vector of solving equation (\in TYPE^number_of_rows)
*  @param acc                  - [in] required accuracy (will not always be achieved)
*  @param max_it               - [in] maximum number of allowed iterations
*
*  @throw exception            - when matrix is not squared
*/
//---------------------------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::iterative_refinement( const input_storage_scheme<TYPE>& ISS,
                                                         TYPE *x,
                                                         TYPE *b,
                                                         double acc,
                                                         size_t max_it
                                                         ) const
{
    if (number_of_columns != number_of_rows)
        throw std::exception("dynamic_storage_scheme<TYPE>::iterative_refinement: matrix is not squared");

    const size_t N = number_of_columns;

    TYPE *d = new TYPE[N];
    TYPE *r = new TYPE[N];
    TYPE *y = new TYPE[N];

    size_t iteration = 0;
    ISS.count_rasidual_vector(x, b, r);
    double v_norm = vector_norm(r, N);
    double new_v_norm;

    // int while condition are contained 2 conditins to stop the calculations,
    // third condition is implemented inside the loop
    // =======================================================================
    while (iteration < max_it && v_norm > acc)
    {
        switch (dynamic_state)
        {
            case LU_DECOMPOSED:
                solve_LU(d, r, y);
                for (size_t i = 0; i < N; ++i)
                    d[i] = x[i] - d[i];
                break;

            case ITERATIVE:
                SOR_iteration(d, b, x);
                break;
        }

        ISS.count_rasidual_vector(d, b, r);
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

    delete[] d;
    delete[] r;
    delete[] y;
}


//--------------------------------------------------------------------------- iterative_preparation
/**
*  Method prepares matrix to SOR iterations
*
*  @throw exception     - when stored matrix has at least one zero diagonal element
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::iterative_preparation(void)
{
    if (dynamic_state != ROL_INIT)
        throw std::exception("dynamic_storage_scheme<TYPE>::iterative_preparation: ROL_INIT is required");

    for (size_t row = 0; row < number_of_rows; ++row)
    {
        const size_t rows_end = HA[row][3];
        for (int idx = HA[row][1]; idx <= (int)rows_end; ++idx)
        {
            const size_t col_number = (size_t)CNLU[idx];
            if (col_number < row && HA[row][2] == idx)
                ++HA[row][2];
            else if (col_number < row)
            {
                TYPE val = ALU[idx];
                CNLU[idx] = CNLU[HA[row][2]];
                ALU[idx] = ALU[HA[row][2]];
                CNLU[HA[row][2]] = col_number;
                ALU[HA[row][2]++] = val;
            }
            else if (col_number == row)
            {
                PIVOT[row] = ALU[idx];
                eliminate_element_in_ROL(idx, row);
                --idx;
            }
        }
        if (PIVOT[row] == 0)
            throw std::exception("iterative_preparation: diagonal element is equal to zero");
    }

    dynamic_state = ITERATIVE;
}

//----------------------------------------------------------------------------------- SOR_iteration
/**
*  Method performs one iteration of SOR method
*
*  @param x                     - [in]  obtained solution
*  @param b                     - [in]  right sight vestor of equation Ax=b
*  @param prev_x                - [out] started approximation
*  @param oemga - (default = 1) - [in]  relaxation parameter 
*                                       default set to 1 (Gauss-Seidel method)
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::SOR_iteration( TYPE *x,
                                                  TYPE *b,
                                                  TYPE *prev_x
                                                  ) const
{
    if (dynamic_state != ITERATIVE)
        std::exception("dynamic_storage_scheme<TYPE>::SOR_iteration: iterative_preparation is needed before");
    if (number_of_columns != number_of_rows)
        std::exception("dynamic_storage_scheme<TYPE>::SOR_iteration: matrix is not squared");

    for (size_t row = 0; row < number_of_rows; ++row)
    {
        x[row] = b[row];
        int idx = HA[row][1];
        for (; idx < HA[row][2]; ++idx)
            x[row] -= (ALU[idx] * x[CNLU[idx]]);
        for (; idx <= HA[row][3]; ++idx)
            x[row] -= (ALU[idx] * prev_x[CNLU[idx]]);
        x[row] *= omega;
        x[row] /= PIVOT[row];
        x[row] += (1 - omega) * prev_x[row];
    }
}

//-------------------------------------------------------------------------- print_sparsity_pattern
/**
*  Method prints sparsity patern of a matrix to the specified file
*
*  @param file_name       - [in] name of a output file
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::print_sparsity_pattern(const char *file_name)
{
std::ofstream outFile;

    // Try to open/create the file
    // ===========================
    try
    {
        outFile.open(file_name);
    }
    catch(std::exception)
    {
        logged_errors += "dynamic_storage_scheme<TYPE>::print_sparsity_pattern: can not open the file";
        logged_errors += std::string(file_name);
        outFile.close();
        return;
    }

    for (size_t row = 0; row < number_of_rows; row++)
    {
        for (size_t col = 0; col < number_of_columns; col++)
        {
            bool found = false;
            for (int idx = HA[HA[row][7]][1]; idx <= HA[HA[row][7]][3]; idx++)
                if (HA[CNLU[idx]][10] == col)
                {
                    outFile << "X";
                    found = true;
                    break;
                }
            if(!found)
                outFile << " ";
        }
        outFile << std::endl;
    }

outFile.close();
}



//***********************************************************************************************//
//                                     PRIVATE MEMBER FUNCTIONS                                  //
//***********************************************************************************************//

//------------------------------------------------------------------------------------ permute_rows
/**
*  Functions permutes rows that are acctualy on pos1 and pos2 positions
*
*  @param pos1          - [in] current position in matrix of first specified row
*  @param pos2          - [in] current position in matrix of second specified row
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void dynamic_storage_scheme<TYPE>::permute_rows( size_t pos1,
                                                        size_t pos2
                                                        )
{
    if (pos1 == pos2)
        return;

    if (pos1 >= number_of_rows || pos2 >= number_of_rows)
        throw std::out_of_range("dynamic_storage_scheme<TYPE>::permute_rows: values out of range");

    HA[HA[pos1][7]][8] = pos2;
    HA[HA[pos2][7]][8] = pos1;

    int val;
    val = HA[pos1][7];
    HA[pos1][7] = HA[pos2][7];
    HA[pos2][7] = val;
}

//--------------------------------------------------------------------------------- permute_columns
/**
*  Functions permutes columns that are acctualy on pos1 and pos2 positions
*
*  @param pos1          - [in] current position in matrix of first specified column
*  @param pos2          - [in] current position in matrix of second specified column
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void dynamic_storage_scheme<TYPE>::permute_columns( size_t pos1,
                                                           size_t pos2
                                                           )
{
    if (pos1 == pos2)
        return;

    if (pos1 >= number_of_columns || pos2 >= number_of_columns)
        throw std::out_of_range("dynamic_storage_scheme<TYPE>::permute_columns: values out of range");

    HA[HA[pos1][9]][10] = pos2;
    HA[HA[pos2][9]][10] = pos1;

    int val;
    val = HA[pos1][9];
    HA[pos1][9] = HA[pos2][9];
    HA[pos2][9] = val;
}

//-------------------------------------------------------------------------------- store_fillin_ROL
/**
*  Function is use for storing new element in row-ordered list
*  function checks many cases in most optimal way
*
*  @param val                              - [in] value of storing element
*  @param row                              - [in] row number of storing element
*  @param col                              - [in] column number of storing element
*  @param garbage_on - (default = true)    - [in] flag indicating if garbage collection
*                                            should be performed
*
*  @return STORING_SUCCESS                 - if element was successfuly stored
*  @return STORING_FAIL                    - if element wasn't successfuly stored
*/
//-------------------------------------------------------------------------------------------------

template <typename TYPE>
int dynamic_storage_scheme<TYPE>::store_fillin_ROL( TYPE val, int origRow,
                                                    int origCol,
                                                    bool garbage_on
                                                    )
{
    // if row isn't empty
    // ==================
    if (HA[origRow][1] != FREE)
    {
        if (CROL < NROL)
        {
            const int after = HA[origRow][3] + 1,
                      before = HA[origRow][1] - 1;
    
            // store element after row package
            // ===============================
            if ((size_t)after < NROL && CNLU[after] == FREE)
            {
                ALU[after] = val;
                CNLU[after] = origCol;
                HA[origRow][3] = after;
                CROL++;
                if ((size_t)after > LROL)
                    LROL = (size_t)after;
            }
            // store element before row package
            // ================================
            else if (before >= 0 && CNLU[before] == FREE)
            {
                const int idx = HA[origRow][2] - 1;
                ALU[before] = ALU[idx];
                CNLU[before] = CNLU[idx];
                ALU[idx] = val;
                CNLU[idx] = origCol;
                HA[origRow][2]--;
                HA[origRow][1] = before;
                CROL++;
            }
            // copy row package to further position
            // ====================================
            else if ((size_t)(HA[origRow][3] - HA[origRow][1] + 2) <= (NROL - LROL - 1))
            {
                int idx;
                const int dist = LROL + 1 - HA[origRow][1];
                for (idx = HA[origRow][1]; idx <= HA[origRow][3]; idx++)
                {
                    ALU[idx + dist] = ALU[idx];
                    CNLU[idx + dist] = CNLU[idx];
                    CNLU[idx] = FREE;
                }
                ALU[idx + dist] = val;
                CNLU[idx + dist] = origCol;
                HA[origRow][1] += dist;
                HA[origRow][2] += dist;
                HA[origRow][3] += (dist + 1);
                CROL++;
                LROL = (size_t)HA[origRow][3];
            }
            // else garbage collection is needed
            // =================================
            // note: in some cases garbage collection cann't help to store element even there is a place to do it,
            // to avoid stack overflow in this cases, this block of code should be executed only once during storing element
            // =============================================================================================================
            else if (garbage_on)
            {
                garbage_collection_in_ROL();
                return store_fillin_ROL(val, origRow, origCol, false);
            }
            else
            {
                logged_errors += "dynamic_storage_scheme<TYPE>::store_fillin_ROL: out of memory\n";
                return STORING_FAIL;
            }
        }
        else
        {
            logged_errors += "dynamic_storage_scheme<TYPE>::store_fillin_ROL: out of memory\n";
            return STORING_FAIL;
        }
    }
    // in case when row is empty
    // =========================
    else
    {
        const size_t idx = LROL + 1;
        if (idx < NROL && CNLU[idx] == FREE)
        {
            HA[origRow][1] = HA[origRow][2] = HA[origRow][3] = idx;
            ALU[idx] = val;
            CNLU[idx] = origCol;
            CROL = 1;
            LROL++;
        }
        else if (garbage_on)
        {
            garbage_collection_in_ROL();
            return store_fillin_ROL(val, origRow, origCol, false);
        }
        else
        {
            logged_errors += "dynamic_storage_scheme<TYPE>::store_fillin_ROL: out of memory\n";
            return STORING_FAIL;
        }
    }

    return STORING_SUCCESS;
}

//----------------------------------------------------------------------  garbage_collection_in_ROL
/**
*  Function performs the organize of the elements in a compact structure
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::garbage_collection_in_ROL(void)
{
    // Mark beginings of row packages by setting numbers less then FREE  in CNLU
    // =========================================================================
    for (size_t row = 0; row < number_of_rows; row++)
    {
        const int k = HA[row][1];
        if (k == FREE)
            continue;
        HA[row][1] = CNLU[k];
        CNLU[k] = - (int)row - 2;
    }
    // Repositioning of all elements
    // =============================
    int NEWPOS = 0;
    for(size_t idx = 0; idx <= LROL; idx++)
    {
         if (CNLU[idx] == FREE)
             continue;
         else if (CNLU[idx] < FREE)
         {
              ALU[NEWPOS] = ALU[idx];
              const int row = - CNLU[idx] - 2;
              CNLU[NEWPOS] = HA[row][1];
              HA[row][1] = NEWPOS;
              HA[row][2] = NEWPOS + HA[row][2] - idx;
              HA[row][3] = NEWPOS + HA[row][3] - idx;
              NEWPOS++;
              continue;
         }
         else
         {
             ALU[NEWPOS] = ALU[idx];
             CNLU[NEWPOS] = CNLU[idx];
             NEWPOS++;
         }
    }
    for (size_t idx = NEWPOS; idx <= LROL; idx++)
        CNLU[idx] = FREE;

    LROL = NEWPOS - 1;
}

//-------------------------------------------------------------------------------- store_fillin_ROL
/**
*  Function is use for storing new element in column-ordered list
*  function checks many cases in most optimal way
*
*  @param row                              - [in] row number of storing element
*  @param col                              - [in] column number of storing element
*  @param garbage_on (default = true)      - [in] flag indicating if garbage collection
*                                            should be performed
*
*  @return STORING_SUCCESS                 - if element was successfuly stored
*  @return STORING_FAIL                    - if element wasn't successfuly stored
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
int dynamic_storage_scheme<TYPE>::store_fillin_COL( int origRow,
                                                    int origCol,
                                                    bool garbage_on
                                                    )
{
    // if column isn't empty
    // =====================
    if (HA[origCol][4] != FREE)
    {
        if (CCOL < NROL)
        {
            const int after = HA[origCol][6] + 1,
                      before = HA[origCol][4] - 1;
    
            // store element after column package
            // ==================================
            if ((size_t)after < NCOL && RNLU[after] == FREE)
            {
                RNLU[after] = origRow;
                HA[origCol][6] = after;
                CCOL++;
                if ((size_t)after > LCOL)
                    LCOL = (size_t)after;
            }
            // store element before column package
            // ===================================
            else if (before >= 0 && RNLU[before] == FREE)
            {
                const int idx = HA[origCol][5] - 1;
                RNLU[before] = RNLU[idx];
                RNLU[idx] = origRow;
                HA[origCol][5]--;
                HA[origCol][4] = before;
                CCOL++;
            }
            // copy column package to further position
            // =======================================
            else if ((size_t)(HA[origCol][6] - HA[origCol][4] + 2) <= (NCOL - LCOL - 1))
            {
                int idx;
                const int dist = LCOL + 1 - HA[origCol][4];
                for (idx = HA[origCol][4]; idx <= HA[origCol][6]; idx++)
                {
                    RNLU[idx + dist] = RNLU[idx];
                    RNLU[idx] = FREE;
                }
                RNLU[idx + dist] = origRow;
                HA[origCol][4] += dist;
                HA[origCol][5] += dist;
                HA[origCol][6] += (dist + 1);
                CCOL++;
                LCOL = (size_t)HA[origCol][6];
            }
            // else garbage collection is needed
            // =================================
            // note: in some cases garbage collection cann't help to store element even there is a place to do it,
            // to avoid stack overflow in this cases, this block of code should be executed only once during storing element
            // =============================================================================================================
            else if (garbage_on)
            {
                garbage_collection_in_COL();
                return store_fillin_COL(origRow, origCol, false);
            }
            else
            {
                logged_errors += "dynamic_storage_scheme<TYPE>::store_fillin_COL: out of memory\n";
                return STORING_FAIL;
            }
        }
        else
        {
            logged_errors += "dynamic_storage_scheme<TYPE>::store_fillin_COL: out of memory\n";
            return STORING_FAIL;
        }
    }
    // in case when column is empty
    // ============================
    else
    {
        const size_t idx = LCOL + 1;
        if (idx < NCOL && RNLU[idx] == FREE)
        {
            HA[origCol][4] = HA[origCol][5] = HA[origCol][6] = idx;
            RNLU[idx] = origRow;
            CCOL = 1;
            LCOL++;
        }
        else if (garbage_on)
        {
            garbage_collection_in_COL();
            return store_fillin_COL(origRow, origCol, false);
        }
        else
        {
            logged_errors += "dynamic_storage_scheme<TYPE>::store_fillin_COL: out of memory\n";
            return STORING_FAIL;
        }
    }

    return STORING_SUCCESS;
}

//----------------------------------------------------------------------  garbage_collection_in_ROL
/**
*  Function performs the organize of the elements in a compact structure
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::garbage_collection_in_COL(void)
{
    // Mark beginings of column packages by setting numbers less then FREE in RNLU
    // ===========================================================================
    for (size_t col = 0; col < number_of_columns; col++)
    {
         const int k = HA[col][4];
         if (k == FREE)
            continue;
         HA[col][4] = RNLU[k];
         RNLU[k] = - (int)col - 2;
    }
    // Repositioning of all elements
    // =============================
    int NEWPOS = 0;
    for(size_t idx = 0; idx <= LCOL; idx++)
    {
         if (RNLU[idx] == FREE)
             continue;
         else if (RNLU[idx] < FREE)
         {
              const int col = - RNLU[idx] - 2;
              RNLU[NEWPOS] = HA[col][4];
              HA[col][4] = NEWPOS;
              HA[col][5] = NEWPOS + HA[col][5] - idx;
              HA[col][6] = NEWPOS + HA[col][6] - idx;
              NEWPOS++;
              continue;
         }
         else
         {
             RNLU[NEWPOS] = RNLU[idx];
             NEWPOS++;
         }
    }
    for (size_t idx = NEWPOS; idx <= LCOL; idx++)
        RNLU[idx] = FREE;

    LCOL = NEWPOS - 1;
}


//----------------------------------------------------------------------- deactivate_element_in_ROL
/**
*  Function bring choosen element to inactive part of specified row package
*
*  @param index                - [in] index in ROL of element we wont to deactive
*  @param row                  - [in] row-package of speccified element
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void dynamic_storage_scheme<TYPE>::deactivate_element_in_ROL( size_t index,
                                                                     size_t row
                                                                     )
{
    TYPE val = ALU[index];
    ALU[index] = ALU[HA[row][2]];
    ALU[HA[row][2]] = val;
    size_t col = CNLU[index];
    CNLU[index] = CNLU[HA[row][2]];
    CNLU[HA[row][2]] = col;
    HA[row][2]++;
}
//------------------------------------------------------------------------ eliminate_element_in_ROL
/**
*  Function deletes choosen element from specified row-package
*
*  @param index                - [in] index in ROL of element we wont to elimate
*  @param row                  - [in] row-package of speccified element
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void dynamic_storage_scheme<TYPE>::eliminate_element_in_ROL( size_t index,
                                                                    size_t row
                                                                    )
{
    ALU[index] = ALU[HA[row][3]];
    CNLU[index] = CNLU[HA[row][3]];
    CNLU[HA[row][3]] = FREE;
    if (LROL == HA[row][3])
        LROL--;
    CROL--;
    HA[row][3]--;
    if (HA[row][3] < HA[row][1])
    {
        HA[row][1] = HA[row][2] = FREE;
        HA[row][3] = EMPTY;
    }
}

//------------------------------------------------------------------------ eliminate_element_in_COL
/**
*  Function deletes choosen element from specified column-package
*
*  @param index                - [in] index in COL of element we wont to elimate
*  @param col                  - [in] column-package of speccified element
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
inline void dynamic_storage_scheme<TYPE>::eliminate_element_in_COL( size_t index,
                                                                    size_t col
                                                                    )
{
    RNLU[index] = RNLU[HA[col][6]];
    RNLU[HA[col][6]] = FREE;
    if (LCOL == HA[col][6])
        LCOL--;
    CCOL--;
    HA[col][6]--;
    if (HA[col][6] < HA[col][4])
    {
        HA[col][4] = HA[col][5] = FREE;
        HA[col][6] = EMPTY;
    }
}

//-------------------------------------------------------------------- delete_column_package_in_COL
/**
*  Function deletes choosen column package in column ordered list
*
*  @param index                - [in] index in ROL of element we wont to elimate
*  @param row                  - [in] row-package of speccified element
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::delete_column_package_in_COL(size_t col)
{
const int col_end = HA[col][6];

    for (int col_idx = HA[col][4]; col_idx <= col_end; col_idx++)
        RNLU[col_idx] = FREE;
    CCOL -= (HA[col][6] - HA[col][4] + 1);
    if (LCOL == HA[col][6])
    {
        if (!CCOL)
            LCOL = FREE;
        else
        {
            for (int idx = HA[col][4]; idx >= 0; idx--)
                if (RNLU[idx] != FREE)
                {
                    LCOL = idx;
                    break;
                }
        }
    }
    HA[col][4] = HA[col][5] = FREE;
    HA[col][6] = EMPTY;
}



//====================================== PIVOTAL STRATEGIES =======================================

//--------------------------------------------------------------- choose_pivot_by_ONE_ROW_SEARCHING
/**
*  Function choose pivot which meets stability condition and its column-package
*  has the smalest number of non-zeroes
*
*  @param stage        - [in] stage of elimination
*  @param mult         - [in] stable parametar
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::choose_pivot_by_ONE_ROW_SEARCHING( size_t stage,
                                                                      double mult
                                                                      )
{
size_t INDEX;
double max_val = 0;
double ABS_VAL = 0;
int MIN_COLCOST = 0;

    const int begin_row = HA[HA[stage][7]][2],
              end_row = HA[HA[stage][7]][3];

    // find the greates in absolute value element in row and store proper data
    // =======================================================================
    for (int idx = begin_row; idx <= end_row; idx++)
    {
        const double abs_val = absolute_value(ALU[idx]);
        if (abs_val > max_val)
        {
            max_val = abs_val;
            INDEX = idx;
            const int col = CNLU[idx];
            MIN_COLCOST = HA[col][6] - HA[col][5];
        }
    }
    if (max_val == 0)
    {
        logged_errors += "dynamic_storage_scheme<TYPE>::choose_pivot_by_ONE_ROW_SEARCHING: some active row contains only zeroes in its active part\n";
        return;
    }

    ABS_VAL = max_val;

    // choose the best element in considering row
    // ==========================================
    for (int idx = begin_row; idx <= end_row; idx++)
    {
        const double abs_val = absolute_value(ALU[idx]);
        const int col = CNLU[idx];
        const int min_colcost = HA[col][6] - HA[col][5];
        if (abs_val * mult >= max_val && (min_colcost < MIN_COLCOST || (min_colcost == MIN_COLCOST && abs_val > ABS_VAL)))
        {
            ABS_VAL = abs_val;
            INDEX = idx;
            MIN_COLCOST = min_colcost;
        }
    }
    // bring element to position (stage, stage)
    // ========================================
    permute_columns(HA[CNLU[INDEX]][10], stage);

    // store pivot and free position in ROL
    // ====================================
    PIVOT[stage] = ALU[INDEX];
    eliminate_element_in_ROL(INDEX, HA[stage][7]);
}

//------------------------------------------------------------------ choose_pivot_by_MARKOWITZ_COST
/**
*  Function choose pivot which markowitz cost is optimal
*
*  @param stage        - [in] stage of elimination
*  @param search       - [in] number of rows to search
*  @param mult         - [in] stable parametar
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::choose_pivot_by_MARKOWITZ_COST( size_t stage,
                                                                   size_t search,
                                                                   double mult
                                                                   )
{
const size_t LastSearch = (stage + search < number_of_rows ? stage + search : number_of_rows);
size_t INDEX;
double ABS_VAL = 0;
int MIN_MCOST = 0;
size_t ROW;

    for (size_t row = stage; row < LastSearch; row++)
    {
        double max_val = 0;
        int index = FREE;

        // get begin and end of the active part of row
        // ===========================================
        const int begin_row = HA[HA[row][7]][2],
                  end_row = HA[HA[row][7]][3];

        // find the greatest absolute value in row
        // =======================================
        for (int idx = begin_row; idx <= end_row; idx++)
        {
            if (idx <= FREE)
            {
                logged_errors += "dynamic_storage_scheme<TYPE>::choose_pivot_by_MARKOWITZ_COST: some active row is empty\n";
                return;
            }
            const double temp_abs_val = absolute_value(ALU[idx]);
            if (temp_abs_val > max_val)
            {
                max_val = temp_abs_val;
                index = idx;
            }
        }
        double abs_val = max_val;
        const int row_mult = end_row - begin_row;
        if (index == FREE)
        {
            logged_errors += "dynamic_storage_scheme<TYPE>::choose_pivot_by_MARKOWITZ_COST: some active row contains only zeroes in its active part\n";
            return;
        }
        int min_mcost = row_mult * (HA[CNLU[index]][6] - HA[CNLU[index]][5]);
        for (int idx = begin_row; idx <= end_row; idx++)
        {
            // if element meets the stability criterion then consider it as a pivot
            // ====================================================================
            const double temp_abs_val = absolute_value(ALU[idx]);
            if (temp_abs_val * mult >= max_val)
            {
                // choose the best element in current row
                // ======================================
                const int mcost = row_mult * (HA[CNLU[idx]][6] - HA[CNLU[idx]][5]);
                if (mcost < min_mcost || (mcost == min_mcost && temp_abs_val > abs_val))
                {
                    abs_val = temp_abs_val;
                    index = idx;
                    min_mcost = mcost;
                }
            }
        }
        // choose the best element in searching area (all considered rows)
        // ===============================================================
        if (min_mcost < MIN_MCOST || (min_mcost == MIN_MCOST && abs_val > ABS_VAL) || row == stage)
        {
            MIN_MCOST = min_mcost;
            ABS_VAL = abs_val;
            ROW = row;
            INDEX = index;
        }
    }
    // bring element to position (stage, stage)
    // ========================================
    permute_rows(ROW, stage);
    permute_columns(HA[CNLU[INDEX]][10], stage);

    // store pivot and free position in ROL
    // ====================================
    PIVOT[stage] = ALU[INDEX];
    eliminate_element_in_ROL(INDEX, HA[stage][7]);

}

//----------------------------------------------------------- choose_pivot_by_FILLIN_MINIMALIZATION
/**
*  Function choose pivot which couses the smalest number of fillins
*
*  @param stage        - [in] stage of elimination
*  @param search       - [in] number of rows to search
*  @param mult         - [in] stable parametar
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::choose_pivot_by_FILLIN_MINIMALIZATION( size_t stage,
                                                                          size_t search,
                                                                          double mult
                                                                          )
{
const size_t LastSearch = (stage + search < number_of_rows ? stage + search : number_of_rows);
size_t INDEX;
double ABS_VAL = 0;
int MIN_FCOST = 0;
size_t ROW;

    for (size_t row = stage; row < LastSearch; row++)
    {
        double max_val = 0;
        int index = FREE;

        // get begin and end of the active part of row
        // ===========================================
        const int begin_row = HA[HA[row][7]][2],
                  end_row = HA[HA[row][7]][3];

        // find the greatest absolute value in row
        // =======================================
        for (int idx = begin_row; idx <= end_row; idx++)
        {
            if (idx <= FREE)
            {
                logged_errors += "dynamic_storage_scheme<TYPE>::choose_pivot_by_FILLIN_MINIMALIZATION: some active row is empty\n";
                return;
            }
            const double temp_abs_val = absolute_value(ALU[idx]);
            if (temp_abs_val > max_val)
            {
                max_val = temp_abs_val;
                index = idx;
            }
        }
        if (index == FREE)
        {
            logged_errors += "dynamic_storage_scheme<TYPE>::choose_pivot_by_FILLIN_MINIMALIZATION: some active row contains only zeroes in its active part\n";
            return;
        }
        int max_index = index;
        double abs_val = max_val;
        int min_fcost = count_fillin_cost(max_index, HA[row][7]);
        for (int idx = begin_row; idx <= end_row; idx++)
        {
            if (idx == max_index)
                continue;

            // if element meets the stability criterion then consider it as a pivot
            // ====================================================================
            const double temp_abs_val = absolute_value(ALU[idx]);
            if (temp_abs_val * mult >= max_val)
            {
                // choose the best element in current row
                // ======================================
                const int fcost = count_fillin_cost(idx, HA[row][7]);
                if (fcost < min_fcost || (fcost == min_fcost && temp_abs_val > abs_val))
                {
                    abs_val = temp_abs_val;
                    index = idx;
                    min_fcost = fcost;
                }
            }
        }
        // choose the best element in searching area (all considered rows)
        // ===============================================================
        if (min_fcost < MIN_FCOST || (min_fcost == MIN_FCOST && abs_val > ABS_VAL) || row == stage)
        {
            MIN_FCOST = min_fcost;
            ABS_VAL = abs_val;
            ROW = row;
            INDEX = index;
        }
    }
    // bring element to position (stage, stage)
    // ========================================
    permute_rows(ROW, stage);
    permute_columns(HA[CNLU[INDEX]][10], stage);

    // store pivot and free position in ROL
    // ====================================
    PIVOT[stage] = ALU[INDEX];
    eliminate_element_in_ROL(INDEX, HA[stage][7]);
}

//--------------------------------------------------------------------------------------- sort_rows
/**
*  Standard quick sort algorithm, used for sorting rows in not decreasing order of
*  non-zero elements
*
*  @param l            - [in] left index of sorting range
*  @param r            - [in] right index of sorting range
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
void dynamic_storage_scheme<TYPE>::sort_rows( int l,
                                              int r
                                              )
{
    if (!r)
        r = number_of_rows - 1;

    const int mid = (l + r) / 2;
    const int v = HA[HA[mid][7]][3] - HA[HA[mid][7]][2];
    int i = l,
        j = r;
    do
    {
        while (i < (int)number_of_rows && HA[HA[i][7]][3] - HA[HA[i][7]][2] < v)
            i++;
        while (j >= 0 && v < HA[HA[j][7]][3] - HA[HA[j][7]][2])
            j--;
        if(i <= j)
        {
            permute_rows(i, j);
            i++;
            j--;
        }
    } while(i < j);
    if(l < j)
        sort_rows(l, j);
    if(i < r)
        sort_rows(i, r);
}


//------------------------------------------------------------------------------- count_fillin_cost
/**
*  Method counts fillin cost on stage "stage" of elimination in assumption
*  that specified element was choosen as a pivot
*
*  @param index            - [in] index of specified element in ROL
*  @param row              - [in] original row number of specified element
*
*  @return                 - fillin cost of element
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
int dynamic_storage_scheme<TYPE>::count_fillin_cost( size_t index,
                                                     size_t row
                                                     )
{
size_t FCOST = 0;
const int col_number = CNLU[index];

    const int row_begin = HA[row][2],
              row_end = HA[row][3];

    // mark proper positions in specjal positions HA[][0]
    // and store potential number of fillins in one eliminated row
    // ===========================================================
    size_t num_elems = 0;
    for (int idx = row_begin; idx <= row_end; ++idx)
        if (idx != index)
        {
            HA[HA[CNLU[idx]][10]][0] = NFREE;
            ++num_elems;
        }
    size_t temp_fcost = num_elems;

    // over all eliminated rows
    // ========================
    const int check_col_end = HA[col_number][6];
    for (int col = HA[col_number][5]; col <= check_col_end; ++col)
    {
        const int row_number = RNLU[col];
        if (row_number != row)
        {
            const int checking_row_end = HA[row_number][3];
            for (int idx = HA[row_number][2]; idx <= checking_row_end; ++idx)
            {
                const int col = CNLU[idx];
                if (col != col_number && HA[HA[col][10]][0] == NFREE)
                    --temp_fcost;
            }
            FCOST += temp_fcost;
            temp_fcost = num_elems;
        }
    }

    // clear marked positions
    // ======================
    for (int idx = row_begin; idx <= row_end; ++idx)
        HA[HA[CNLU[idx]][10]][0] = FREE;

    return FCOST;
}
//---------------------------------------- FRIEND FUNCTIONS ---------------------------------------

//-------------------------------------------------------------------------------------- operator<<
/**
*  Standard outstream operator
*
*  @param out                   - [in/out] standard outstream object
*  @param DSS                   - [in] dynamic storage scheme
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE2>
std::ostream& operator<<( std::ostream& out,
                          const dynamic_storage_scheme<TYPE2>& DSS
                          )
{
    const size_t manip_int = 3;
    const size_t manip_idx = 4;
    const size_t manip_double = 10;

    out << "\\=================== DYNCAMIC STORAGE SCHEME ===================" << std::endl
        << "number of rows       :" << DSS.number_of_rows << std::endl
        << "number of columns    :" << DSS.number_of_columns << std::endl
        << "dynamic state        :";
    switch (DSS.dynamic_state)
        {
        case ROL_INIT: out << "ROL_INIT\n\n"; break;
        case COL_INIT: out << "COL_INIT\n\n"; break;
        case ITERATIVE: out << "ITERATIVE\n\n"; break;
        case LU_DECOMPOSED: out << "LU_DECOMPOSED\n\n"; break;
        case QR_DECOMPOSED: out << "QR_DECOMPOSED\n\n"; break;
        }

    out << "PIVOT: [PIVOT[idx],idx]" << std::endl;
    for (size_t i = 0; i < DSS.NHA; i++)
        out << "[" << DSS.PIVOT[i] << "," << i << "]";
    out << std::endl << std::endl;
    out << "\\----------------- ROW ORDERED LIST -----------------" << std::endl
        << "size of row-ordered list                  :" << DSS.NROL << std::endl
        << "last not free position in ROL             :" << DSS.LROL << std::endl
        << "number of non-zeros actualy stored in ROL :" << DSS.CROL << std::endl << std::endl;

    // Print memory in usege
    // =====================
    out << "memory usage   :";
    for (size_t Pos = 0; Pos < DSS.NROL; Pos++)
    {
        if (DSS.CNLU[Pos] != FREE)
            out << "o";
        else 
            out << ".";
    }
    out << std::endl;

    // Print indexes
    // =============
    out << "indexes % 1    :";
    for (size_t Pos = 0; Pos < DSS.NROL; Pos++)
        out << Pos % 10;
    out << std::endl;
    out << "indexes % 10   :";
    for (size_t Pos = 0; Pos < DSS.NROL; Pos++)
    {
        if (Pos % 10 == 0)
            out << (Pos / 10) % 10;
        else
            out << " ";
    }
    out << std::endl;

    // Print part of integrity table responsible for ROL
    // =================================================
    for (size_t Index = 1; Index <= 3; Index++)
    {
        // Print ROL markers
        // =================
        out << "H[.][" << Index << "] % 10   :";
        for (size_t Pos = 0; Pos < DSS.NROL; Pos++)
        {
            bool SignPrinted = false;
            for (size_t RowNum = 0; RowNum < DSS.number_of_rows; RowNum++)
            {
                if (DSS.HA[RowNum][Index] == Pos &&
                    !(Index == 2 && DSS.HA[RowNum][2] > DSS.HA[RowNum][3] && (size_t)DSS.HA[RowNum][2] < DSS.NROL && DSS.CNLU[DSS.HA[RowNum][2]] != FREE))
                {
                    out << (RowNum % 10);
                    SignPrinted = true;
                }
            }
            if (!SignPrinted)
                out << " ";
        }
        out << std::endl;
    }
    out << std::endl;
    // Print markers
    // =============
    out << "row =   ";
    for (size_t row = 0; row < DSS.number_of_rows; row++)
        out << std::setw(manip_idx) << row;
    out << std::endl;

    // Print part of integrity table responsible for ROL
    // =================================================
    for (size_t Index = 1; Index <= 3; Index++)
    {
        // Print ROL markers
        // =================
        out << "H[.][" << Index << "]:";
        for (size_t row = 0; row < DSS.number_of_rows; row++)
            out << std::setw(manip_idx) << DSS.HA[row][Index];
        out << std::endl;
    }

    // Print stored in ROL values, their column numbers and positions in ROL
    // =====================================================================
    out << std::endl << "(original row number, its current position):";
    if (DSS.dynamic_state != COL_INIT)
        out << "[ALU,CNLU,idx],..." << std::endl;
    else
        out << "[CNLU,idx],..." << std::endl;
    for (size_t row = 0; row < DSS.number_of_rows; row++)
    {
        out << "(" << DSS.HA[row][7] << "," << row << "): ";
        if (DSS.HA[DSS.HA[row][7]][1] == FREE || DSS.HA[DSS.HA[row][7]][2] == FREE || DSS.HA[DSS.HA[row][7]][3] == EMPTY)
            continue;
        for (int idx = DSS.HA[DSS.HA[row][7]][1]; idx <= DSS.HA[DSS.HA[row][7]][3]; idx++)
        {
            if (DSS.dynamic_state != COL_INIT)
                out << "[" << std::setw(manip_double) << DSS.ALU[idx] <<","<< std::setw(manip_int) << DSS.CNLU[idx] << "," << std::setw(manip_int) << idx << "]";
            else
                out << "[" << std::setw(manip_int) << DSS.CNLU[idx] << "," << std::setw(manip_int) << idx << "]";
        }
        out << std::endl;
    }

    // Print permutations of rows
    // ==========================
    out << std::endl << "row permutations:" << std::endl;
    out << "         ";
    for (size_t row = 0; row < DSS.number_of_rows; row++)
        out << std::setw(manip_int) << row;
    out << std::endl << "HA[][7]: ";
    for (size_t row = 0; row < DSS.number_of_rows; row++)
        out << std::setw(manip_int) << DSS.HA[row][7];
    out << std::endl << "HA[][8]: ";
    for (size_t row = 0; row < DSS.number_of_rows; row++)
        out << std::setw(manip_int) << DSS.HA[row][8];

    out << std::endl << std::endl;
    out << "\\---------------- COLUMN ORDERED LIST ----------------" << std::endl
        << "size of column-ordered list               :" << DSS.NCOL << std::endl
        << "last not free position in COL             :" << DSS.LCOL << std::endl
        << "number of non-zeros actualy stored in COL :" << DSS.CCOL << std::endl << std::endl;

    // Print memory in usege
    // =====================
    out << "memory usage   :";
    for (size_t Pos = 0; Pos < DSS.NCOL; Pos++)
    {
        if (DSS.RNLU[Pos] != FREE)
            out << "o";
        else 
            out << ".";
    }
    out << std::endl;

    // Print indexes
    // =============
    out << "indexes % 1    :";
    for (size_t Pos = 0; Pos < DSS.NCOL; Pos++)
        out << Pos % 10;
    out << std::endl;
    out << "indexes % 10   :";
    for (size_t Pos = 0; Pos < DSS.NCOL; Pos++)
    {
        if (Pos % 10 == 0)
            out << (Pos / 10) % 10;
        else
            out << " ";
    }
    out << std::endl;

    // Print part of integrity table responsible for COL
    // =================================================
    for (size_t Index = 4; Index <= 6; Index++)
    {
        // Print COL markers
        // =================
        out << "H[.][" << Index << "] % 10   :";
        for (size_t Pos = 0; Pos < DSS.NCOL; Pos++)
        {
            bool SignPrinted = false;
            for (size_t ColNum = 0; ColNum < DSS.number_of_columns; ColNum++)
            {
                if (DSS.HA[ColNum][Index] == Pos &&
                    !(Index == 5 && DSS.HA[ColNum][5] > DSS.HA[ColNum][6] && (size_t)DSS.HA[ColNum][5] < DSS.NCOL && DSS.RNLU[DSS.HA[ColNum][5]] != FREE))
                {
                    out << (ColNum % 10);
                    SignPrinted = true;
                }
            }
            if (!SignPrinted)
                out << " ";
        }
        out << std::endl;
    }
    out << std::endl;
    // Print markers
    // =============
    out << "col =   ";
    for (size_t col = 0; col < DSS.number_of_columns; col++)
        out << std::setw(manip_idx) << col;
    out << std::endl;

    // Print part of integrity table responsible for COL
    // =================================================
    for (size_t Index = 4; Index <= 6; Index++)
    {
        // Print ROL markers
        // =================
        out << "H[.][" << Index << "]:";
        for (size_t col = 0; col < DSS.number_of_columns; col++)
            out << std::setw(manip_idx) << DSS.HA[col][Index];
        out << std::endl;
    }

    // Print stored in COL values and positions in COL
    // ===============================================
    out << std::endl << "(original column number, its current position):";
    if (DSS.dynamic_state != COL_INIT)
        out << "[RNLU,idx],..." << std::endl;
    else
        out << "[ALU,RNLU,idx],..." << std::endl;
    for (size_t col = 0; col < DSS.number_of_columns; col++)
    {
        if (DSS.dynamic_state != COL_INIT)
            out << "(" << std::setw(manip_int) << DSS.HA[col][9] << "," << std::setw(manip_int) << col << ")";
        else
        {
            out << std::setw(manip_double + 1) << " " << "(" << std::setw(manip_int) << DSS.HA[col][9] << "," << std::setw(manip_int) << col << ")";
        }
    }
    out << std::endl;

    size_t line = 0;
    bool print = true;
    while (print && line < DSS.number_of_rows)
    {
        print = false;
        for (size_t col = 0; col < DSS.number_of_columns; col++)
        {
            if (DSS.HA[DSS.HA[col][9]][4] == FREE)
            {
                out << std::setw(2 * manip_int + 3) << " ";
                continue;
            }
            int idx = DSS.HA[DSS.HA[col][9]][4] + line;
            if (idx <= DSS.HA[DSS.HA[col][9]][6])
            {
                print = true;
                if (DSS.dynamic_state != COL_INIT)
                    out << "[" << std::setw(manip_int) << DSS.RNLU[idx] << "," << std::setw(manip_int) << idx << "]";
                else
                    out << "[" <<  std::setw(manip_double) << DSS.ALU[idx] << "," << std::setw(manip_int) << DSS.RNLU[idx] << "," << std::setw(manip_int) << idx << "]";
            }
            else
            {
                if (DSS.dynamic_state != COL_INIT)
                    out << std::setw(2 * manip_int + 3) << " ";
                else
                    out << std::setw(manip_double + 2 * manip_int + 4) << " ";
            }
        }
        out << std::endl;
        line++;
    }

    // Print permutations of rows
    // ==========================
    out << std::endl << "column permutations:" << std::endl;
    out << "         ";
    for (size_t col = 0; col < DSS.number_of_columns; col++)
        out << std::setw(manip_int) << col;
    out << std::endl << "HA[][9]: ";
    for (size_t col = 0; col < DSS.number_of_columns; col++)
        out << std::setw(manip_int) << DSS.HA[col][9];
    out << std::endl << "HA[][10]:";
    for (size_t col = 0; col < DSS.number_of_columns; col++)
        out << std::setw(manip_int) << DSS.HA[col][10];

    // Print working parts of integrity table
    // ======================================
    out << std::endl << std::endl << "\\---------------- WORKING PART OF INTEGRITY TABLE  ----------------" << std::endl;
    out << "         ";
    for (size_t col = 0; col < DSS.number_of_columns; col++)
        out << std::setw(manip_int) << col;
    out << std::endl << "HA[][0]: ";
    for (size_t col = 0; col < DSS.NHA; col++)
        out << std::setw(manip_int) << DSS.HA[col][0];

    out << std::endl << std::endl << "\\============================ MATRIX RECONSTRUCTION ============================";
    if (DSS.dynamic_state != COL_INIT)
    {
        for (size_t row = 0; row < DSS.number_of_rows; row++)
        {
            out << std::endl;
            for (size_t col = 0; col < DSS.number_of_columns; col++)
            {
                bool printed = false;
                for (int idx = DSS.HA[DSS.HA[row][7]][1]; idx <= DSS.HA[DSS.HA[row][7]][3]; idx++)
                {
                    if (idx <= FREE)
                        break;
                    if (DSS.CNLU[idx] == DSS.HA[col][9])
                    {
                        out << std::setw(manip_double) << DSS.ALU[idx];
                        printed = true;
                        break;
                    }
                }
                if (!printed && row == col && absolute_value(DSS.PIVOT[row]) != 0)
                    out << std::setw(manip_double) << DSS.PIVOT[row];
                else if (!printed)
                    out << std::setw(manip_double) << 0;
            }
        }
    }
    else
    {
    }

    return out;
}



/*************************************************************************************************/
/*                                                                                               */
/*                                    TESTING FUNCTIONS                                          */
/*                                                                                               */
/*************************************************************************************************/
//---------------------------------------------------------------------------- check_integrity_test
/**
*  Function checks the consistency of the data contained in dynamic_storage_scheme
*/
//-------------------------------------------------------------------------------------------------
template <typename TYPE>
int dynamic_storage_scheme<TYPE>::check_integrity_test(void) const
{
int a, b, c;
int ret_val = 0;

    cout << "\n=== Checking dynamic_storage_scheme ===\n";

    //------------------- ROW ORDERED LIST CHECKING -------------------
    // row ordered list checking
    // =========================
    a = b = c = 0;
    for (size_t idx = 0; idx < NROL; idx++)
    {
        if (CNLU[idx] == FREE)
            a++;
        else if (CNLU[idx] < 0)
            b++;
        else
            c++;
    }
    if (a + c != NROL)
    {
        cout << "memory error in ROL: " << b << " pbits" << std::endl;
        ret_val = 2;
    }

    // check more detailed stored data
    // ===============================
    for (size_t row = 0; row < number_of_rows; row++)
    {
        if (HA[row][1] > HA[row][2] || HA[row][1] > HA[row][3])
            if(HA[row][1] >= 0 && HA[row][1] < (int)NROL && CNLU[HA[row][1]] != FREE)
            {
                cout << "ROL integrity error: row " << row << " should be empty or begin set to FREE" << std::endl;
                ret_val = 2;
            }

        for (int idx = HA[row][1]; idx <= HA[row][3]; idx++)
        {
            if (idx <= FREE)
            {
                cout << "ROL warning: row " << row << " is empty" << std::endl;
                ret_val = 1;
                if (HA[row][1] != FREE || HA[row][2] != FREE || HA[row][3] != EMPTY)
                {
                    cout << "ROL error: 1, 2 and 3 pointers of row " << row << " should be set to FREE or EMPTY" << std::endl;
                    ret_val = 2;
                    break;
                }
            }
            else if (CNLU[idx] < 0 || ALU[idx] == 0)
            {
                cout << "ROL integrity lack error in row " << row << std::endl;
                ret_val = 2;
            }
        }
    }

    //------------------ COLUMN ORDERED LIST CHECKING -----------------
    // column ordered list checking
    // ============================
    a = b = c = 0;
    for (size_t idx = 0; idx < NCOL; idx++)
    {
        if (RNLU[idx] == FREE)
            a++;
        else if (RNLU[idx] < 0)
            b++;
        else
            c++;
    }
    if (a + c != NCOL)
    {
        cout << "memory error in COL: " << b << " pbits" << std::endl;
        ret_val = 2;
    }

    // check more detailed stored data
    // ===============================
    for (size_t col = 0; col < number_of_columns; col++)
    {
        if (HA[col][4] > HA[col][5] || HA[col][4] > HA[col][6])
            if(HA[col][4] != FREE)
            {
                cout << "COL integrity error: col " << col << " should be empty or begin set to FREE" << std::endl;
                ret_val = 2;
            }

        for (int idx = HA[col][4]; idx <= HA[col][6]; idx++)
        {
            if (idx <= FREE)
            {
                cout << "COL warning: column " << col << " is empty" << std::endl;
                ret_val = 1;
                if (HA[col][4] != FREE || HA[col][5] != FREE || HA[col][6] != EMPTY)
                {
                    cout << "COL error: 4, 5 and 6 pointers of column " << col << " should be set to FREE or EMPTY" << std::endl;
                    ret_val = 2;
                    break;
                }
            }
            if (RNLU[idx] < 0)
            {
                cout << "COL integrity lack error in col " << col << std::endl;
                ret_val = 2;
            }
        }
    }

    // check permutations parts
    // ========================
    for (size_t row = 0; row < number_of_rows; row++)
        if (HA[HA[row][7]][8] != row || HA[HA[row][8]][7] != row)
        {
            cout << "permutations error row: " << row << std::endl;
            ret_val = 2;
        }
    for (size_t col = 0; col < number_of_columns; col++)
        if (HA[HA[col][9]][10] != col || HA[HA[col][10]][9] != col)
        {
            cout << "permutations error row: " << col << std::endl;
            ret_val = 2;
        }

    // check the membership every non-zero elem in ROL to any row-package
    // ==================================================================
    for (size_t idx = 0; idx < NROL; idx++)
    {
        if (CNLU[idx] > FREE)
        {
            bool checked = false;
            for (size_t row = 0; row < number_of_rows; row++)
                if (HA[row][1] <= (int)idx && (int)idx <= HA[row][3])
                {
                    checked = true;
                    break;
                }
            if (!checked)
            {
                cout << "ROL error: element: " << idx << " doesn't belong to any row package" << std::endl;
                ret_val = 2;
            }
        }
        else if (CNLU[idx] == FREE)
        {
            bool checked = true;
            size_t row;
            for (row = 0; row < number_of_rows; row++)
                if (HA[row][1] <= (int)idx && (int)idx <= HA[row][3])
                {
                    checked = false;
                    break;
                }
            if (!checked)
            {
                cout << "ROL error: empty element: " << idx << " belongs to row: " << row << std::endl;
                ret_val = 2;
            }
        }
        else
        {
            cout << "ROL error: invalid value in CNLU: index = " << idx << std::endl;
            ret_val = 2;
        }
    }

    // check the membership every non-zero elem in COL to any column-package
    // =====================================================================
    for (size_t idx = 0; idx < NCOL; idx++)
    {
        if (RNLU[idx] > FREE)
        {
            bool checked = false;
            for (size_t col = 0; col < number_of_columns; col++)
                if (HA[col][4] <= (int)idx && (int)idx <= HA[col][6])
                {
                    checked = true;
                    break;
                }
            if (!checked)
            {
                cout << "COL error: element: " << idx << " doesn't belong to any column package" << std::endl;
                ret_val = 2;
            }
        }
        else if (RNLU[idx] == FREE)
        {
            bool checked = true;
            size_t col;
            for (col = 0; col < number_of_columns; col++)
                if (HA[col][4] <= (int)idx && (int)idx <= HA[col][6])
                {
                    checked = false;
                    break;
                }
            if (!checked)
            {
                cout << "COL error: empty element: " << idx << " belongs to column: " << col << std::endl;
                ret_val = 2;
            }
        }
        else
        {
            cout << "COL error: invalid value in RNLU: index = " << idx << std::endl;
            ret_val = 2;
        }
    }

    int *check_tab = NULL;
    try
    {
        check_tab = new int[NHA];
    }
    catch(std::bad_alloc)
    {
        delete[] check_tab;
        cout << "warning: some data can not be checked becouse of bad_alloc:check_tab exception" << std::endl;
        ret_val = 1;
        cout << "=== dynamic_storage_scheme checked ===" << std::endl;
        return ret_val;
    }
    // check duplication elements in row packages in ROL
    // =================================================
    for (size_t row = 0; row < number_of_rows; row++)
    {
        for (size_t col = 0; col < number_of_columns; col++)
            check_tab[col] = col;
        for (int idx = HA[row][1]; idx <= HA[row][3]; idx++)
        {
            if (idx <= FREE) // in case when row is empyt
                break;
            const int col = CNLU[idx];
            if (col >= 0 && col < (int)number_of_columns && check_tab[col] == col)
                check_tab[col] = FREE;
            else if (col >= 0 && col < (int)number_of_columns && check_tab[col] == FREE)
            {
                cout << "ROL error: duplicate element: column: " << col << " in original row: " << row << std::endl;
                ret_val = 2;
            }
            else
            {
                cout << "ROL error: range error ocured or wrong column number in row: " << row << std::endl;
                ret_val = 2;
            }
        }
    }
    // check duplication elements in column packages in COL
    // ====================================================
    for (size_t col = 0; col < number_of_columns; col++)
    {
        for (size_t row = 0; row < number_of_rows; row++)
            check_tab[row] = row;
        for (int idx = HA[col][4]; idx <= HA[col][6]; idx++)
        {
            if (idx <= FREE) // in case when col is empyt
                break;
            const int row = RNLU[idx];
            if (row >= 0 && row < (int)number_of_rows && check_tab[row] == row)
                check_tab[row] = FREE;
            else if (row >= 0 && row < (int)number_of_columns && check_tab[row] == FREE)
            {
                cout << "COL error: duplicate element: row: " << row << " in original column: " << col << std::endl;
                ret_val = 2;
            }
            else
            {
                cout << "ROL error: range error ocured or wrong row number in column: " << col << std::endl;
                ret_val = 2;
            }
        }
    }

    delete[] check_tab;

    cout << "=== dynamic_storage_scheme checked ===" << std::endl;

    return ret_val;
}

#endif


