#include <filesystem>
#include <gtest/gtest.h>

#include <sparse_matrices.hpp>
#include <matrix_rand_generators.hpp>


using namespace std;

constexpr double eps_float{ 1e-3 };
constexpr double eps_double{ 1e-10 };

static size_t g_test_id{ 0 };


template< typename T >
class non_singular_linear_equation : public ::testing::Test
{
protected:
	input_storage_scheme< T > ISS;
	unique_ptr< dynamic_storage_scheme< T > > DSS;

	vector< T > b, x, r;
	T low_val{ 0.01 }, high_val{ 10.0 };

	filesystem::path test_dir{ "test_prints" };

	virtual input_storage_scheme< T > generate_matrix() {	return generate_ISS< T >( get_mx_size(), get_mx_size(), true, get_zero_proportion(), low_val, high_val );	}
	virtual size_t get_mx_size() { return 200; }
	virtual size_t get_zero_proportion() { return 40; }
	virtual DYNAMIC_STATE get_init_type() { return DYNAMIC_STATE::ROL_INIT; };

	void SetUp() override
	{
		++g_test_id;

		if( !filesystem::exists( test_dir ) )
			filesystem::create_directory( test_dir );

		string s_pattern_file = test_dir.string() + "/sparsity_pattern_test_" + to_string( g_test_id ) + ".txt";

		if( std::is_same_v< real_type< T >::type, double> )
		{
			low_val = 0.00001;
			high_val = 10000.0;
		}

		// generate random non singular matrix
		ISS = generate_matrix();
		// permute order in input schemy as all algorithms should works independetly for order of adding
		EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );
		// create dynamic scheme
		EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< T > >( ISS, 2, 0.8, get_init_type() ) );
		// just print sparsity pattern
		EXPECT_NO_THROW( DSS->print_sparsity_pattern( s_pattern_file.c_str() ) );

		// allocation of vectors of equation Ax = b
		// ========================================
		b.resize( get_mx_size(), 0.0 );
		x.resize( get_mx_size(), 0.0 );
		r.resize( get_mx_size(), 0.0 ); // residual vector ( inaccuracy vector )

		random_vector_values( &b, low_val, high_val );
	}

	void TearDown() override
	{
		// function that checks integrity of dynamic scheme
		EXPECT_EQ( DSS->check_integrity_test(), 0 );
		// amount of non zeros after dynamic operations can be seen here
		size_t non_zeros{ 0 };
		EXPECT_NO_THROW( non_zeros = DSS->get_non_zeros_amount() );
		// imporve the result by iterative refinemnt
		EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );
		// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
		EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

		// check inacurracy no some accuarcy level
		auto acuracy{ l2_norm( r ) / l2_norm( b ) };
		if( std::is_same_v< real_type< T >::type, float> )
			EXPECT_LE( acuracy, eps_float );
		if( std::is_same_v< real_type< T >::type, double> )
			EXPECT_LE( acuracy, eps_double );
	}
};

using test_types = ::testing::Types< float, double, complex< float >, complex< double > >;

TYPED_TEST_SUITE( non_singular_linear_equation, test_types );

TYPED_TEST( non_singular_linear_equation, LU_decomposition_MarkowitzCost_NONE )
{
	// decompose A=LU using Gauss elimination
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, get_mx_size(), 1.0, 0.0, LD_PREPARATION::NONE ) );
}

TYPED_TEST( non_singular_linear_equation, LU_decomposition_MarkowitzCost_SORT )
{	
	// decompose A=LU using Gauss elimination
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, 10, 1.0, 0.0, LD_PREPARATION::SORT ) );
}

TYPED_TEST( non_singular_linear_equation, LU_decomposition_MarkowitzCost_AMD )
{
	// decompose A=LU using Gauss elimination
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, 10, 1.0, 0.0, LD_PREPARATION::AMD ) );
}

TYPED_TEST( non_singular_linear_equation, LU_decomposition_FillinMin_AMD )
{
	// decompose A=LU using Gauss elimination
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, 1, 1.0, 0.0, LD_PREPARATION::AMD ) );
}





template< typename T >
class non_singular_linear_equation_COL_INIT : public non_singular_linear_equation< T >
{
protected:
	size_t get_mx_size() override { return 8; }
	size_t get_zero_proportion() override { return 3; }
	DYNAMIC_STATE get_init_type() override { return DYNAMIC_STATE::COL_INIT; };
	void TearDown() override {};
};

TYPED_TEST_SUITE( non_singular_linear_equation_COL_INIT, test_types );

TYPED_TEST( non_singular_linear_equation_COL_INIT, QR_decomposition_sort_cols )
{
	// decompose A=QR using Householder alghoritm
	EXPECT_NO_THROW( DSS->QR_decomposition( LD_PREPARATION::AMD ) );

}





template< typename T >
class non_singular_linear_equation_strong_diag : public non_singular_linear_equation< T >
{
protected:
	input_storage_scheme< T > generate_matrix() override
	{
		return generate_ISS_with_strong_diagonal< T >( get_mx_size(), get_mx_size(), get_zero_proportion(), low_val, high_val );
	}
};

TYPED_TEST_SUITE( non_singular_linear_equation_strong_diag, test_types );

TYPED_TEST( non_singular_linear_equation_strong_diag, SOR_iterations )
{
	EXPECT_NO_THROW( DSS->iterative_preparation() );
}


template< typename T >
class small_matrix_just_for_print_ROL : public non_singular_linear_equation< T >
{
protected:
	virtual size_t get_mx_size() override { return 10; }
	virtual size_t get_zero_proportion() override { return 2; }
	virtual void TearDown() override {}
};

TYPED_TEST_SUITE( small_matrix_just_for_print_ROL, test_types );

TYPED_TEST( small_matrix_just_for_print_ROL, print_matrix )
{
	string s_pattern_file = test_dir.string() + "/print_scheme_ROL_test_" + to_string( g_test_id ) + ".txt";
	DSS->print_scheme_to_file( s_pattern_file.c_str() );
}


template< typename T >
class small_matrix_just_for_print_COL : public small_matrix_just_for_print_ROL< T >
{
protected:
	DYNAMIC_STATE get_init_type() override { return DYNAMIC_STATE::COL_INIT; };
};

TYPED_TEST_SUITE( small_matrix_just_for_print_COL, test_types );

TYPED_TEST( small_matrix_just_for_print_COL, print_matrix )
{
	string s_pattern_file = test_dir.string() + "/print_scheme_COL_test_" + to_string( g_test_id ) + ".txt";
	DSS->print_scheme_to_file( s_pattern_file.c_str() );
}
