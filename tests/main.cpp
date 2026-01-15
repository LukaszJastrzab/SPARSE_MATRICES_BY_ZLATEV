#include <gtest/gtest.h>

#include <sparse_matrices.hpp>
#include <matrix_rand_generators.hpp>


using namespace std;


TEST( NonSingularLinearEquation, LU_decomposition )
{
	const size_t MxSize = 100;

	auto ISS = GenerateISS<double>( MxSize, MxSize, true, 1, 0.00001, 100000.0 );

	unique_ptr< dynamic_storage_scheme<double> > DSS;
	EXPECT_NO_THROW( DSS = make_unique< dynamic_storage_scheme<double> >( ISS, 10, 0.8 ) );
	EXPECT_NO_THROW( DSS->LU_decomposition( PIVOTAL_STRATEGY::MARKOWITZ_COST, 8, 5, numeric_limits<double>::min(), true ) );

	// allocation of vectors
	// =====================
	vector<double> b( MxSize, 0.0 );
	vector<double> x( MxSize, 0.0 );

	RandomVectorValues( &b, 0.00001, 10000.0 );

	// solve equation to obtain first aproximation of the solution
	// ===========================================================
	EXPECT_NO_THROW( DSS->solve_LU( x.data(), b.data() ) );
	EXPECT_NO_THROW( DSS->iterative_refinement( ISS, x.data(), b.data(), 0.0000000000000000001, 20 ) );


	//EXPECT_EQ( 2 + 2, 4 );
}