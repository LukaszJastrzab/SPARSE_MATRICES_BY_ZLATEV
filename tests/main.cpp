#include <gtest/gtest.h>

#include <sparse_matrices.hpp>
#include <matrix_rand_generators.hpp>


using namespace std;


TEST( non_singular_linear_equation_real, LU_decomposition_markowitz )
{
	const size_t mx_size = 50;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 1, 0.00001, 100000.0 );

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
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, mx_size, 1.0, numeric_limits< double >::min(), true ) );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// but in some very hard example this condition will not be met
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( non_singular_linear_equation_complex, LU_decomposition_markowitz )
{
	const size_t mx_size = 50;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< complex< double > >( mx_size, mx_size, true, 1, 0.00001, 100000.0 );

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
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, mx_size, 1.0, numeric_limits< double >::min(), true ) );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// but in some very hard example this condition will not be met
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( non_singular_linear_equation_real, LU_decomposition_fillin_min )
{
	const size_t mx_size = 50;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 1, 0.00001, 100000.0 );

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
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< double >::min(), true ) );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// but in some very hard example this condition will not be met
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( non_singular_linear_equation_complex, LU_decomposition_fillin_min )
{
	const size_t mx_size = 50;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< complex< double > >( mx_size, mx_size, true, 1, 0.00001, 100000.0 );

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
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::FILLIN_MINIMALIZATION, mx_size, 1.0, numeric_limits< double >::min(), true ) );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// but in some very hard example this condition will not be met
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( non_singular_linear_equation_real, LU_decomposition_one_row_search )
{
	const size_t mx_size = 50;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< double >( mx_size, mx_size, true, 1, 0.00001, 100000.0 );

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
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::ONE_ROW_SEARCHING, mx_size, 1.0, numeric_limits< double >::min(), true ) );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// but in some very hard example this condition will not be met
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( non_singular_linear_equation_complex, LU_decomposition_one_row_search )
{
	const size_t mx_size = 50;

	// Input Store Scheme
	// ==================
	auto ISS = generate_ISS< complex< double > >( mx_size, mx_size, true, 1, 0.00001, 100000.0 );

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
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::ONE_ROW_SEARCHING, mx_size, 1.0, numeric_limits< double >::min(), true ) );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// but in some very hard example this condition will not be met
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( SOR_linear_equation_real, iterative_preparation )
{
	const size_t mx_size = 50;

	// as SOR method is not alwas convergent
	// thus we will create matrix with strong diagonal domination in row
	// |a_ii| = sum_j (|a_ij|)
	// this will ensure SOR method convergence
	// =================================================================
	auto ISS = generate_ISS_with_strong_diagonal(mx_size, mx_size, 1, 0.0001, 10000.0 );

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

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( SOR_linear_equation_complex, iterative_preparation )
{
	const size_t mx_size = 50;

	// as SOR method is not alwas convergent
	// thus we will create matrix with strong diagonal domination in row
	// |a_ii| = sum_j (|a_ij|)
	// this will ensure SOR method convergence
	// =================================================================
	auto ISS = generate_ISS_with_strong_diagonal< complex< double > >( mx_size, mx_size, 1, 0.0001, 10000.0 );

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

	// test if iterative_preparation keeps the rules of handling with dynamic_scheme
	// just for test purposes
	//==========================================
	EXPECT_EQ( DSS->check_integrity_test(), 0 );

	// allocation of vectors of equation Ax = b
	// ========================================
	vector< complex< double > > b( mx_size, complex < double >( 0.0 ) );
	vector< complex< double > > x( mx_size, complex < double >( 0.0 ) );
	vector< complex< double > > r( mx_size, complex < double >( 0.0 ) ); // residual vector ( inaccuracy vector )

	random_vector_values( &b, complex< double >( 0.00001 ), complex< double >( 10000.0 ) );

	// imporve the result by iterative refinemnt
	// =========================================
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x, b, numeric_limits< double >::min(), 1000 ) );

	// calculate residual vector, i.e. r = b - Ax ( in "on-paper" solution r should be zeros )
	// =======================================================================================
	EXPECT_NO_THROW( ISS.count_rasidual_vector( x, b, r ) );

	// calculate residual vector norm - to see the inaccuracy
	// ======================================================
	double norm = l2_norm( r );

	// check inacurracy no some accuarcy level
	// due to not big matrix we should gave even more precise result
	// =============================================================
	EXPECT_TRUE( norm >= 0 && norm < 0.00000001 );
}

TEST( PrintScheme, real_matrix )
{
	const size_t mx_size = 20;

	auto ISS = generate_ISS_with_strong_diagonal( mx_size, mx_size, 1, 0.0001, 10000.0 );

	unique_ptr< dynamic_storage_scheme< double > > DSS;

	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme< double > >( ISS, 1, 1 ) );

	std::ofstream file( "PrintScheme_real_matrix.txt" );

	EXPECT_TRUE( file );

	EXPECT_NO_THROW( file << *DSS );

	file.close();
}