#ifndef MATRIX_SCHEMES_RANDOM_GENERATORS
#define MATRIX_SCHEMES_RANDOM_GENERATORS

#include <algorithm>
#include <random>
#include <vector>

#include <sparse_matrices.hpp>

std::random_device rd;
std::mt19937 g( rd() );

std::uniform_int_distribution< int > bool_rng( 0, 1 );

template < typename T >
void random_vector_values( std::vector< T >* vec, T min_val, T max_val )
{
	if( min_val <= 0 || max_val <= 0 )
		throw;
	
	std::uniform_real_distribution< double > elem_rng( min_val, max_val );

	for( auto& item : *vec )
		item = ( bool_rng( g ) ? -1.0 : 1.0 )* elem_rng( g );
}

template <>
void random_vector_values( std::vector< std::complex< double> >* vec, std::complex< double> min_val, std::complex< double> max_val )
{
	if( min_val.real() <= 0 || max_val.real() <= 0 )
		throw;

	std::uniform_real_distribution< double > elem_rng( min_val.real(), max_val.real() );

	for( auto& item : *vec )
	{
		double sign_real = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );
		double sign_img = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );	

		item = std::complex< double >( sign_real * elem_rng( g ), sign_img * elem_rng( g ) );
	}
}

template < typename T >
input_storage_scheme< T > generate_ISS( size_t mx_rows, size_t mx_cols, bool non_singularity, int zeros_proportion, T min_val, T max_val )
{
	if( min_val <= 0 || max_val <= 0 || zeros_proportion < 0 )
		throw;

	size_t mx_min_size = std::min( mx_rows, mx_cols );

	input_storage_scheme< double > ISS( mx_rows, mx_cols );

	std::vector< size_t > singular_mx_permuts;
	singular_mx_permuts.reserve( mx_min_size );
	for( size_t idx{ 0 }; idx < mx_min_size; ++idx )
		singular_mx_permuts.push_back( idx );

	std::shuffle( singular_mx_permuts.begin(), singular_mx_permuts.end(), g );

	// randomizers
	// ===========
	std::uniform_int_distribution< int > int_rng( 0, zeros_proportion );
	std::uniform_real_distribution< double > elem_rng( min_val, max_val );

	for( size_t row{ 0 }; row < mx_rows; ++row )
		for( size_t col{ 0 }; col < mx_cols; ++col )

			if( ( int_rng( g ) && ( non_singularity || singular_mx_permuts[ row ] != col ) ) ||
				( non_singularity && singular_mx_permuts[ row ] == col ) )
			{
				T sign = static_cast< T >( bool_rng( g ) ? -1.0 : 1.0 );
				EXPECT_NO_THROW( ISS.add_element( sign * elem_rng( g ), row, col ) );
			}

	return ISS;
}

template <>
input_storage_scheme< std::complex< double > > generate_ISS( size_t mx_rows,
															size_t mx_cols,
															bool non_singularity,
															int zeros_proportion,
															std::complex< double > min_val,
															std::complex< double > max_val )
{
	if( min_val.real() <= 0 || max_val.real() <= 0 || zeros_proportion < 0 )
		throw;

	size_t mx_min_size = std::min( mx_rows, mx_cols );

	input_storage_scheme< std::complex< double > > ISS( mx_rows, mx_cols );

	std::vector< size_t > singular_mx_permuts;
	singular_mx_permuts.reserve( mx_min_size );
	for( size_t idx{ 0 }; idx < mx_min_size; ++idx )
		singular_mx_permuts.push_back( idx );

	std::shuffle( singular_mx_permuts.begin(), singular_mx_permuts.end(), g );

	// randomizers
	// ===========
	std::uniform_int_distribution< int > int_rng( 0, zeros_proportion );
	std::uniform_real_distribution< double > elem_rng( min_val.real(), max_val.real() );

	for( size_t row{ 0 }; row < mx_rows; ++row )
		for( size_t col{ 0 }; col < mx_cols; ++col )

			if( ( int_rng( g ) && ( non_singularity || singular_mx_permuts[ row ] != col ) ) ||
				( non_singularity && singular_mx_permuts[ row ] == col ) )
			{
				double sign_real = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );
				double sign_img  = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );

				std::complex< double > val( sign_real * elem_rng( g ), sign_img * elem_rng( g ) );

				EXPECT_NO_THROW( ISS.add_element( val, row, col ) );
			}

	return ISS;
}

template < typename T >
input_storage_scheme< T > generate_ISS_with_strong_diagonal( size_t mx_rows, size_t mx_cols, int zeros_proportion, T min_val, T max_val )
{
	if( min_val <= 0 || max_val <= 0 || zeros_proportion < 0 )
		throw;

	input_storage_scheme< double > ISS( mx_rows, mx_cols );

	// randomizers
	// ===========
	std::uniform_int_distribution< int > int_rng( 0, zeros_proportion );
	std::uniform_real_distribution< double > elem_rng( min_val, max_val );

	for( size_t col{ 0 }; col < mx_cols; ++col )
	{
		T abs_sum{};

		for( size_t row{ 0 }; row < mx_rows; ++row )
			if( row != col && int_rng( g ) )
			{
				T sign = static_cast< T >( bool_rng( g ) ? -1.0 : 1.0 );
				T val = elem_rng( g );
				abs_sum += val;

				EXPECT_NO_THROW( ISS.add_element( sign * val, row, col ) );
			}

		// ensure strong diagonal
		// ======================
		EXPECT_NO_THROW( ISS.add_element( abs_sum + 0.01, col, col ) );
	}

	return ISS;
}

template <>
input_storage_scheme< std::complex< double > > generate_ISS_with_strong_diagonal( size_t mx_rows,
																				  size_t mx_cols,
																				  int zeros_proportion,
																				  std::complex< double > min_val,
																				  std::complex< double > max_val )
{
	if( min_val.real() <= 0 || max_val.real() <= 0 || zeros_proportion < 0 )
		throw;

	input_storage_scheme< std::complex< double > > ISS( mx_rows, mx_cols );

	// randomizers
	// ===========
	std::uniform_int_distribution< int > int_rng( 0, zeros_proportion );
	std::uniform_real_distribution< double > elem_rng( min_val.real(), max_val.real() );

	for( size_t col{ 0 }; col < mx_cols; ++col )
	{
		double abs_sum{ 0 };

		for( size_t row{ 0 }; row < mx_rows; ++row )
			if( row != col && int_rng( g ) )
			{
				double sign_real = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );
				double sign_img = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );

				std::complex< double > val( sign_real * elem_rng( g ), sign_img * elem_rng( g ) );
				abs_sum += std::abs( val );

				EXPECT_NO_THROW( ISS.add_element( val, row, col ) );
			}

		// ensure strong diagonal
		// ======================
		double sign_real = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );
		double sign_img = static_cast< double >( bool_rng( g ) ? -1.0 : 1.0 );

		abs_sum += 0.1;
		auto real_img_val = sqrt( abs_sum * abs_sum / 2 );
		EXPECT_NO_THROW( ISS.add_element( std::complex< double >( sign_real * real_img_val, sign_img * real_img_val ), col, col ) );
	}

	return ISS;
}

#endif