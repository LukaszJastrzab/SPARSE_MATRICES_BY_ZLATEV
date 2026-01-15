#ifndef MATRIX_SCHEMES_RANDOM_GENERATORS
#define MATRIX_SCHEMES_RANDOM_GENERATORS

#include <algorithm>
#include <random>
#include <vector>

#include <sparse_matrices.hpp>

std::random_device rd;
std::mt19937 g( rd() );

std::uniform_int_distribution<int> boolRng( 0, 1 );

template <typename T>
void RandomVectorValues( std::vector<T>* vec, T minVal, T maxVal )
{
	if( minVal <= 0 || maxVal <= 0 )
		throw;
	
	std::uniform_real_distribution<double> elemRng( minVal, maxVal );

	for( auto& item : *vec )
		item = ( boolRng( g ) ? -1.0 : 1.0 )* elemRng( g );
}


template <typename T>
input_storage_scheme<T> GenerateISS( size_t MxRows, size_t MxCols, bool nonSingularity, int zerosProportion, T minVal, T maxVal )
{
	if( minVal <= 0 || maxVal <= 0 || zerosProportion < 0 )
		throw;

	size_t MxMinSize = min( MxRows, MxCols );

	input_storage_scheme<double> ISS( MxRows, MxCols );

	std::vector<size_t> singularMxPermuts;
	singularMxPermuts.reserve( MxMinSize );
	for( size_t idx{ 0 }; idx < MxMinSize; ++idx )
		singularMxPermuts.push_back( idx );

	std::shuffle( singularMxPermuts.begin(), singularMxPermuts.end(), g );

	// randomizers
	// ===========
	std::uniform_int_distribution<int> intRng( 0, zerosProportion );
	std::uniform_real_distribution<double> elemRng( minVal, maxVal );

	for( size_t row{ 0 }; row < MxRows; ++row )
		for( size_t col{ 0 }; col < MxCols; ++col )

			if( ( intRng( g ) && ( nonSingularity || singularMxPermuts[ row ] != col ) ) ||
				( nonSingularity && singularMxPermuts[ row ] == col ) )
			{
				T sign = static_cast<T>( boolRng( g ) ? -1.0 : 1.0 );
				EXPECT_NO_THROW( ISS.add_element( sign * elemRng( g ), row, col ) );
			}

	return ISS;
}

template <typename T>
input_storage_scheme<T> GenerateISSWithStrongDiagonal( size_t MxRows, size_t MxCols, int zerosProportion, T minVal, T maxVal )
{
	if( minVal <= 0 || maxVal <= 0 || zerosProportion < 0 )
		throw;

	size_t MxMinSize = min( MxRows, MxCols );

	input_storage_scheme<double> ISS( MxRows, MxCols );

	// randomizers
	// ===========
	std::uniform_int_distribution<int> intRng( 0, zerosProportion );
	std::uniform_real_distribution<double> elemRng( minVal, maxVal );

	for( size_t col{ 0 }; col < MxCols; ++col )
	{
		T AbsSum{};

		for( size_t row{ 0 }; row < MxRows; ++row )
			if( row != col && intRng( g ) )
			{
				T sign = static_cast< T >( boolRng( g ) ? -1.0 : 1.0 );
				T val = elemRng( g );
				AbsSum += val;

				EXPECT_NO_THROW( ISS.add_element( sign * val, row, col ) );
			}

		// ensure strong diagonal
		// ======================
		EXPECT_NO_THROW( ISS.add_element( AbsSum + 0.01, col, col ) );
	}

	return ISS;
}

#endif