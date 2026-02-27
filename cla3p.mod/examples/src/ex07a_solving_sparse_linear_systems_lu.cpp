/**
 * @example ex07a_solving_sparse_linear_systems_lu.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/sparse.hpp>
#include <cla3p/linsol.hpp>
#include <cla3p/algebra.hpp>

/*-----------------------------------------------------*/
static cla3p::csr::RdMatrix DefaultSparseMatrix()
{
	static cla3p::int_t  rowptr[] = {0, 3, 5, 8, 11, 13};
	static cla3p::int_t  colidx[] = {0,  1,  3,  0, 1, 2, 3, 4,  0, 2, 3, 1,  4};
	static cla3p::real_t values[] = {1, -1, -3, -2, 5, 4, 6, 4, -4, 2, 7, 8, -5};

	return cla3p::csr::RdMatrix(5, 5, rowptr, colidx, values, false);
}
/*-----------------------------------------------------*/
int main()
{
	const cla3p::csr::RdMatrix A = DefaultSparseMatrix();
	const cla3p::dns::RdVector b = cla3p::dns::RdVector::random(5);
	const cla3p::dns::RdMatrix B = cla3p::dns::RdMatrix::random(5, 3);
	cla3p::dns::RdVector x; // x will be created in solve
	cla3p::dns::RdMatrix X(5, 3); // Preallocate space for X

	cla3p::PardisoLU<cla3p::csr::RdMatrix> luSolver;

	/*
	 * Perform analysis & symbolic decomposition on A.
	 */
	luSolver.analysis(A);

	/*
	 * Decompose A into LU product.
	 */
	luSolver.decompose(A);

	{
		/*
		 * Single column (vector) rhs.
		 * Calculate x = (A^{-1} * b).
		 */
		luSolver.solve(b, x);
		std::cout << "Dense Vector rhs::Absolute Error: "
			      << (b - A * x).evaluate().normOne() << std::endl;
	}

	{
		/*
		 * Multiple column (matrix) rhs.
		 * Calculate X = (A^{-1} * B).
		 */
		luSolver.solve(B, X);
		std::cout << "Dense Matrix rhs::Absolute Error: "
			      << (B - A * X).evaluate().normOne() << std::endl;
	}

	return 0;
}
