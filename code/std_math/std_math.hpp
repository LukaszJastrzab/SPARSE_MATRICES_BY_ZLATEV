#ifndef STANDARD_MATH
#define STANDARD_MATH

#include <iostream>
#include <fstream>

namespace std_math {

//----------------------------------------------------------------------- STANDARD IN/OUT FUNCTIONS
/// Function to saving vector in file
//===================================
template <typename TYPE>
void save_vector( TYPE *v,
                  size_t n,
                  const char* file_name
                  )
{
std::ofstream out_file;

    // Try to open/create the file
    // ===========================
    try
    {
        out_file.open(file_name);
    }
    catch(std::exception)
    {
        out_file.close();
        return;
    }
    out_file << n;
    out_file << std::endl;

    for (size_t i = 0; i < n; ++i)
    {
        out_file << v[i];
        out_file << std::endl;
    }

out_file.close();
}

/// Function to loading vector from file
//======================================
template <typename TYPE>
void load_vector( TYPE *v,
                  size_t n,
                  const char* file_name
                  )
{
std::ifstream input_file;
char buffer[128];

    // Open input file
    // ===============
    input_file.open(file_name);

    if (input_file.is_open())
    {
        for (size_t i = 0; i < n; ++i)
        {
            // =============================================== //
            // ADD YOUR OWN CODE TO STORING SPECIFIED OBJECTS  //
            // =============================================== //
        }
        input_file.close();

        // =======================================================
        throw std::exception("std_math::load_vector: fail: not allowed specialized function");
        // becouse the specialization should be implemented in this particular moment
    }
    else
        throw std::exception("std_math::load_vector: can not open the file");
}
/// specialisation for double
template <>
void load_vector( double *v,
                  size_t n,
                  const char* file_name
                  )
{
std::ifstream input_file;
char buffer[128];

    // Open input file
    // ===============
    input_file.open(file_name);

    if (input_file.is_open())
    {
        size_t i = 0;
        while (!input_file.eof() && i < n)
        {
            input_file >> buffer;
            v[i++] = atof(buffer);
        }
        input_file.close();
    }
    else
        throw std::exception("std_math::load_vector: can not open the file");
}

// class of complex number
// =======================
template <typename TYPE>
class complex_number
{
public:
    TYPE real, imag;

    complex_number(TYPE _real, TYPE _imag)
        :
        real(_real),
        imag(_imag)
    {
    }
    complex_number(TYPE _real = 0)
        :
        real(_real),
        imag(0)
    {
    }

    double abs(void)
    {
        return sqrt(static_cast<double>(real * real + imag * imag));
    }

    template <typename TYPE>
    bool operator== (const TYPE& _real)
    {
        if (imag != 0 || real != _real)
            return false;
        else return true;
    }
    template <typename TYPE>
    bool operator!= (TYPE _real)
    {
        if (imag != 0 || _real != real)
            return true;
        else return false;
    }
    template <typename TYPE>
    void operator= (const TYPE& _real)
    {
        real = _real;
        imag = 0;
    }

private:

    template <typename TYPE>
    friend complex_number<TYPE> operator-(complex_number<TYPE> x);

    template <typename TYPE>
    friend complex_number<TYPE> operator+ (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2);
    template <typename TYPE>
    friend complex_number<TYPE> operator- (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2);
    template <typename TYPE>
    friend complex_number<TYPE> operator* (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2);
    template <typename TYPE>
    friend complex_number<TYPE> operator/ (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2);
};

template <typename TYPE>
complex_number<TYPE> operator-(complex_number<TYPE> x)
{
    complex_number<TYPE> val(-x.real,-x.imag);
    return val;
}

template <typename TYPE>
complex_number<TYPE> operator+ (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2)
{
    complex_number<TYPE> val(l1.real + l2.real, l1.imag + l2.imag);
    return val;
}
template <typename TYPE>
complex_number<TYPE> operator- (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2)
{
    complex_number<TYPE> val(l1.real - l2.real, l1.imag - l2.imag);
    return val;
}

template <typename TYPE>
complex_number<TYPE> operator* (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2)
{
    complex_number<TYPE> val(l1.real * l2.real - l1.imag * l2.imag , l1.real * l2.imag + l1.imag * l2.real);
    return val;
}

template <typename TYPE>
complex_number<TYPE> operator/ (const complex_number<TYPE>& l1, const complex_number<TYPE>& l2)
{
    TYPE div = l2.real * l2.real + l2.imag * l2.imag;
    complex_number<TYPE> val((l1.real * l2.real + l1.imag * l2.imag) / div, (l1.imag * l2.real - l1.real * l2.imag) / div);

    return val;
}




// template definition of the absolute value for difeerent numbers
// ===============================================================
template <typename TYPE>
double absolute_value(const TYPE& x)
{
    return x >= 0 ? x : -x;
}
template <>
double absolute_value(const complex_number<double>& x)
{
    return sqrt(x.real * x.real + x.imag * x.imag);
}

};
#endif