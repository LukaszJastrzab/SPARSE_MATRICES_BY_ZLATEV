#include <gtest/gtest.h>

#include <sparse_matrices.hpp>
#include <matrix_rand_generators.hpp>


using namespace std;

constexpr double eps_float{ 1e-3 };
constexpr double eps_double{ 1e-10 };


TEST( non_singular_linear_equation_real_float, LU_decomposition_markowitz )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< float >( mx_size, mx_size, true, 15, 0.01f, 100.0f );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< float > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< float > >( ISS, 10, 0.8 ) );

	// It prints sparcity pattern
	// ==========================
	EXPECT_NO_THROW( DSS->print_sparsity_pattern( "test1_sparsity_pattern.txt" ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, mx_size, 1.0, numeric_limits< float >::min(), LD_PREPARATION::SORT ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< float > b( mx_size, 0.0 );
	vector< float > x( mx_size, 0.0 );
	vector< float > r( mx_size, 0.0 ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, 0.01f, 100.0f );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_float );
}

TEST( non_singular_linear_equation_real_double, LU_decomposition_markowitz )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 15, 0.00001, 100000.0 );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< double > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< double > >( ISS, 10, 0.8 ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, mx_size, 1.0, numeric_limits< double >::min(), LD_PREPARATION::SORT ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< double > b( mx_size, 0.0 );
	vector< double > x( mx_size, 0.0 );
	vector< double > r( mx_size, 0.0 ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, 0.00001, 10000.0 );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
}

TEST( non_singular_linear_equation_complex_float, LU_decomposition_markowitz )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< complex< float > >( mx_size, mx_size, true, 15, 0.01f, 100.0f );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< complex< float > > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< float > > >( ISS, 10, 0.8 ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, mx_size, 1.0, numeric_limits< float >::min(), LD_PREPARATION::NONE ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< complex< float > > b( mx_size, complex < float >( 0.0 ) );
	vector< complex< float > > x( mx_size, complex < float >( 0.0 ) );
	vector< complex< float > > r( mx_size, complex < float >( 0.0 ) ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, complex < float >( 0.01 ), complex < float >( 100.0 ) );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_float );
}

TEST( non_singular_linear_equation_complex_double, LU_decomposition_markowitz )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< complex< double > >( mx_size, mx_size, true, 15, 0.00001, 100000.0 );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< complex< double > > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< double > > >( ISS, 10, 0.8 ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, mx_size, 1.0, numeric_limits< double >::min(), LD_PREPARATION::NONE ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< complex< double > > b( mx_size, complex < double >( 0.0 ) );
	vector< complex< double > > x( mx_size, complex < double >( 0.0 ) );
	vector< complex< double > > r( mx_size, complex < double >( 0.0 ) ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, complex < double >( 0.00001 ), complex < double >( 10000.0 ) );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
}

TEST( non_singular_linear_equation_real_float, LU_decomposition_fillin_min )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< float >( mx_size, mx_size, true, 15, 0.01f, 100.0f );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< float > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< float > >( ISS, 10, 0.8 ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< float >::min(), LD_PREPARATION::AMD ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< float > b( mx_size, 0.0 );
	vector< float > x( mx_size, 0.0 );
	vector< float > r( mx_size, 0.0 ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, 0.01f, 100.0f );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_float );
}

TEST( non_singular_linear_equation_real_double, LU_decomposition_fillin_min )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 15, 0.00001, 100000.0 );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< double > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< double > >( ISS, 10, 0.8 ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< double >::min(), LD_PREPARATION::NONE ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< double > b( mx_size, 0.0 );
	vector< double > x( mx_size, 0.0 );
	vector< double > r( mx_size, 0.0 ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, 0.00001, 10000.0 );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
}

TEST( non_singular_linear_equation_complex_float, LU_decomposition_fillin_min )
{
	const size_t mx_size = 100;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< complex< float > >( mx_size, mx_size, true, 15, 0.01, 100.0 );

	// This function is redundant in computation
	// here we wont to verify if order of adding elements to ISS make a difference
	// !!!ORDER OF ADDING ELEMENTS TO ISS SHOULD NOT MAKE A DIFFERENCE FOR CALULACTIONS!!!
	// ===================================================================================
	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< complex< float > > > DSS;

	// Create dynamic scheme with some additional memmory for dinamic operations
	// =========================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< float > > >( ISS, 10, 0.8 ) );

	// Decompose matric to triangular (lower and uuper) factor L and U
	// this is done using Gauss elimination with picotal strategy
	// ===============================================================
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< float >::min(), LD_PREPARATION::SORT ) );

	// test if LU_decomposition keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< complex< float > > b( mx_size, complex < float >( 0.0 ) );
	vector< complex< float > > x( mx_size, complex < float >( 0.0 ) );
	vector< complex< float > > r( mx_size, complex < float >( 0.0 ) ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, complex < float >( 0.01 ), complex < float >( 100.0 ) );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x, b ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_float );
}

TEST( non_singular_linear_equation_complex_double, LU_decomposition_fillin_min )
{
	const size_t mx_size = 20;

	for( int i{ 0 }; i < 50; ++i )
	{
		auto ISS = generate_ISS< complex< double > >( mx_size, mx_size, true, 10, 0.00001, 100000.0 );

		EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

		unique_ptr< dynamic_storage_scheme< complex< double > > > DSS;

		EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< double > > >( ISS, 1, 1.8 ) );

		EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< double >::min(), LD_PREPARATION::AMD ) );

		EXPECT_EQ( DSS->check_integrity_test(), 0 );

		vector< complex< double > > b( mx_size, complex < double >( 0.0 ) );
		vector< complex< double > > x( mx_size, complex < double >( 0.0 ) );
		vector< complex< double > > r( mx_size, complex < double >( 0.0 ) ); // residual vector ( inaccuracy vector )

		random_vector_values( &b, complex < double >( 0.00001 ), complex < double >( 10000.0 ) );

		// solve equation to obtain first aproximation of the solution
		// ===========================================================
		EXPECT_NO_THROW( DSS->solve_LU( x, b ) );
		EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );
		EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

		EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
	}
}

TEST( SOR_linear_equation_real_float, iterative_preparation )
{
	const size_t mx_size = 100;

	// as SOR method is not alwas convergent
	// thus we will create matrix with strong diagonal domination in row
	// |a_ii| = sum_j (|a_ij|)
	// this will ensure SOR method convergence
	// =================================================================
	auto ISS = generate_ISS_with_strong_diagonal( mx_size, mx_size, 15, 0.01f, 100.0f );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< float > > DSS;

	// Create dynamic scheme with no additional memorry as iterative method will be used
	// =================================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< float > >( ISS, 1, 1 ) );

	// Set some relaxation parameter from range (0,2) 
	// ==============================================
	EXPECT_NO_THROW( DSS->set_relaxation_parameter( 1.0 ) );

	// make some setups in dynamic stroge scheme that allows to make SOR iteration efficient
	// =====================================================================================
	EXPECT_NO_THROW( DSS->iterative_preparation() );

	// test if iterative_preparation keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< float > b( mx_size, 0.0 );
	vector< float > x( mx_size, 0.0 );
	vector< float > r( mx_size, 0.0 ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, 0.01f, 100.0f );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_float );
}

TEST( SOR_linear_equation_real_double, iterative_preparation )
{
	const size_t mx_size = 100;

	// as SOR method is not alwas convergent
	// thus we will create matrix with strong diagonal domination in row
	// |a_ii| = sum_j (|a_ij|)
	// this will ensure SOR method convergence
	// =================================================================
	auto ISS = generate_ISS_with_strong_diagonal(mx_size, mx_size, 15, 0.0001, 10000.0 );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< double > > DSS;

	// Create dynamic scheme with no additional memorry as iterative method will be used
	// =================================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< double > >( ISS, 1, 1 ) );

	// Set some relaxation parameter from range (0,2) 
	// ==============================================
	EXPECT_NO_THROW( DSS->set_relaxation_parameter( 1.0 ) );

	// make some setups in dynamic stroge scheme that allows to make SOR iteration efficient
	// =====================================================================================
	EXPECT_NO_THROW( DSS->iterative_preparation() );

	// test if iterative_preparation keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< double > b( mx_size, 0.0 );
	vector< double > x( mx_size, 0.0 );
	vector< double > r( mx_size, 0.0 ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, 0.00001, 10000.0 );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
}

TEST( SOR_linear_equation_complex_float, iterative_preparation )
{
	const size_t mx_size = 100;

	// as SOR method is not alwas convergent
	// thus we will create matrix with strong diagonal domination in row
	// |a_ii| = sum_j (|a_ij|)
	// this will ensure SOR method convergence
	// =================================================================
	auto ISS = generate_ISS_with_strong_diagonal< complex< float > >( mx_size, mx_size, 15, 0.01f, 100.0f );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< complex< float > > > DSS;

	// Create dynamic scheme with no additional memorry as iterative method will be used
	// =================================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< float > > >( ISS, 1, 1 ) );

	// Set some relaxation parameter from range (0,2) 
	// ==============================================
	EXPECT_NO_THROW( DSS->set_relaxation_parameter( 1.0 ) );

	// make some setups in dynamic stroge scheme that allows to make SOR iteration efficient
	// =====================================================================================
	EXPECT_NO_THROW( DSS->iterative_preparation() );

	// test if iterative_preparation keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< complex< float > > b( mx_size, complex < float >( 0.0 ) );
	vector< complex< float > > x( mx_size, complex < float >( 0.0 ) );
	vector< complex< float > > r( mx_size, complex < float >( 0.0 ) ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, complex< float >( 0.01 ), complex< float >( 100.0 ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< float >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// check inacurracy no some accuarcy level
	// =======================================
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_float );
}

TEST( SOR_linear_equation_complex_double, iterative_preparation )
{
	const size_t mx_size = 100;

	// as SOR method is not alwas convergent
	// thus we will create matrix with strong diagonal domination in row
	// |a_ii| = sum_j (|a_ij|)
	// this will ensure SOR method convergence
	// =================================================================
	auto ISS = generate_ISS_with_strong_diagonal< complex< double > >( mx_size, mx_size, 15, 0.0001, 10000.0 );

	// Dynamic Storage Scheme
	// ======================
	unique_ptr< dynamic_storage_scheme< complex< double > > > DSS;

	// Create dynamic scheme with no additional memorry as iterative method will be used
	// =================================================================================
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< double > > >( ISS, 1, 1 ) );

	// Set some relaxation parameter from range (0,2) 
	// ==============================================
	EXPECT_NO_THROW( DSS->set_relaxation_parameter( 1.0 ) );

	// make some setups in dynamic stroge scheme that allows to make SOR iteration efficient
	// =====================================================================================
	EXPECT_NO_THROW( DSS->iterative_preparation() );

	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	vector< complex< double > > b( mx_size, complex < double >( 0.0 ) );
	vector< complex< double > > x( mx_size, complex < double >( 0.0 ) );
	vector< complex< double > > r( mx_size, complex < double >( 0.0 ) );

	random_vector_values( &b, complex< double >( 0.00001 ), complex< double >( 10000.0 ) );

	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );

	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
}

TEST( PrintScheme_double, real_matrix )
{
	const size_t mx_size = 20;
	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 15, 0.00001, 100000.0 );
	unique_ptr< dynamic_storage_scheme< double > > DSS;
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< double > >( ISS, 1, 1 ) );
	std::ofstream file( "PrintScheme_real_matrix.txt" );
	EXPECT_TRUE( file );
	EXPECT_NO_THROW( file << *DSS );
	file.close();
}

TEST( PrintScheme_double, complex_matrix )
{
	const size_t mx_size = 20;
	auto ISS = generate_ISS< complex< double > >( mx_size, mx_size, true, 10, 0.00001, 100000.0 );
	unique_ptr< dynamic_storage_scheme< complex< double > > > DSS;
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< complex< double > > >( ISS, 1, 1 ) );
	std::ofstream file( "PrintScheme_complex_matrix.txt" );
	EXPECT_TRUE( file );
	EXPECT_NO_THROW( file << *DSS );
	file.close();
}


TEST( non_singular_linear_equation_real_double, LU_decomposition_fillin_minimalization_compare )
{
	const size_t mx_size = 100;

	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 3, 0.00001, 100000.0 );

	EXPECT_NO_THROW( permute_input_matrix_elements_test( &ISS ) );

	unique_ptr< dynamic_storage_scheme< double > > DSS1, DSS2;

	EXPECT_NO_THROW( DSS1 = make_unique< dynamic_storage_scheme< double > >( ISS, 3, 0.8 ) );
	EXPECT_NO_THROW( DSS2 = make_unique< dynamic_storage_scheme< double > >( ISS, 3, 0.8 ) );

	EXPECT_NO_THROW( DSS1->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< double >::min(), LD_PREPARATION::SORT ) );
	EXPECT_NO_THROW( DSS2->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< double >::min(), LD_PREPARATION::AMD ) );

	EXPECT_EQ( DSS1->check_integrity_test(), 0 );
	EXPECT_EQ( DSS2->check_integrity_test(), 0 );

	auto c1 = DSS1->get_non_zeros_amount();
	auto c2 = DSS2->get_non_zeros_amount();

	vector< double > b( mx_size, 0.0 );
	vector< double > x( mx_size, 0.0 );
	vector< double > r( mx_size, 0.0 );

	random_vector_values( &b, 0.00001, 10000.0 );

	EXPECT_NO_THROW( DSS1->solve_LU( x, b ) );
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );

	EXPECT_NO_THROW( DSS2->solve_LU( x, b ) );
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );
	EXPECT_LE( l2_norm( r ) / l2_norm( b ), eps_double );
}
