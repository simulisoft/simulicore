/**
 * @example ex07c_solving_sparse_linear_systems_ldlt.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/sparse.hpp>
#include <cla3p/linsol.hpp>
#include <cla3p/algebra.hpp>

/*-----------------------------------------------------*/
static cla3p::csr::RdMatrix DefaultSymmetricSparseMatrix()
{
	static cla3p::int_t  rowptr[] = {0, 3, 5, 7,  8, 9};
	static cla3p::int_t  colidx[] = {0,  1,  3, 1, 4, 2, 3, 3,  4};
	static cla3p::real_t values[] = {1, -2, -4, 5, 8, 4, 2, 7, -5};

	return cla3p::csr::RdMatrix(5, 5, rowptr, colidx, values, false, cla3p::Property::SymmetricUpper());
}
/*-----------------------------------------------------*/
int main()
{
	const cla3p::csr::RdMatrix A = DefaultSymmetricSparseMatrix();
	const cla3p::dns::RdVector b = cla3p::dns::RdVector::random(5);
	const cla3p::dns::RdMatrix B = cla3p::dns::RdMatrix::random(5, 3);
	cla3p::dns::RdVector x; // x will be created in solve
	cla3p::dns::RdMatrix X(5, 3); // Preallocate space for X

	cla3p::PardisoLDLt<cla3p::csr::RdMatrix> ldltSolver;

	/*
	 * Perform analysis & symbolic decomposition on A.
	 */
	ldltSolver.analysis(A);

	/*
	 * Decompose A into LDL' product.
	 */
	ldltSolver.decompose(A);

	{
		/*
		 * Single column (vector) rhs.
		 * Calculate x = (A^{-1} * b).
		 */
		ldltSolver.solve(b, x);
		std::cout << "Dense Vector rhs::Absolute Error: "
			      << (b - A * x).evaluate().normOne() << std::endl;
	}

	{
		/*
		 * Multiple column (matrix) rhs.
		 * Calculate X = (A^{-1} * B).

		 */
		ldltSolver.solve(B, X);
		std::cout << "Dense Matrix rhs::Absolute Error: "
			      << (B - A * X).evaluate().normOne() << std::endl;
	}

	return 0;
}
