#include <gtest/gtest.h>
#include <algorithm>
#include <random>

#include <sparse_matrices.hpp>


using namespace std;

random_device rd;
mt19937 g(rd());


template <typename T>
input_storage_scheme< T > GenerateISS(size_t MxRows, size_t MxCols, bool nonSingularity, int zerosProportion, T minVal, T maxVal)
{
    if (minVal <= 0 || maxVal <= 0 || zerosProportion < 0)
        throw;

    size_t MxMinSize = min( MxRows, MxCols );

    input_storage_scheme< double > ISS( MxRows, MxCols );

    vector< size_t > singularMxPermuts;
    singularMxPermuts.reserve( MxMinSize );
    for ( size_t idx{ 0 }; idx < MxMinSize; ++idx )
        singularMxPermuts.push_back( idx );

    shuffle( singularMxPermuts.begin(), singularMxPermuts.end(), g );

    uniform_int_distribution< int > boolRng( 0, nonSingularity );
    uniform_real_distribution< double > elemRng( minVal, maxVal );

    for ( size_t row{ 0 }; row < MxRows; ++row )
        for ( size_t col{ 0 }; col < MxCols; ++col )

            if ( ( boolRng( g ) && ( nonSingularity || singularMxPermuts[row] != col ) ) ||
                ( nonSingularity && singularMxPermuts[row] == col ) )
            {
                double sign = (boolRng(g) ? -1.0 : 1.0);
                EXPECT_NO_THROW( ISS.add_element(sign * elemRng(g), row, col) );
            }

    return ISS;
}


TEST( SampleTest, AlwaysPasse )
{
    const size_t MxSize = 100;

    auto ISS = GenerateISS< double >( MxSize, MxSize, true, 1, 0.00001, 100000.0 );

    unique_ptr< dynamic_storage_scheme< double > > DSS;
    EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< double > >( ISS, 10, 0.8 ) );
    EXPECT_NO_THROW( DSS->LU_decomposition( MARKOWITZ_COST, 8, 5, numeric_limits<double>::min(), true ) );

    //EXPECT_EQ( 2 + 2, 4 );
}