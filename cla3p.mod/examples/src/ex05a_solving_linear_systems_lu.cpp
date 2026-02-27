/**
 * @example ex05a_solving_linear_systems_lu.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/linsol.hpp>
#include <cla3p/algebra.hpp>

int main()
{
	const cla3p::dns::RdMatrix A = cla3p::dns::RdMatrix::random(5, 5);
	const cla3p::dns::RdVector b = cla3p::dns::RdVector::random(5);
	const cla3p::dns::RdMatrix B = cla3p::dns::RdMatrix::random(5, 3);

	cla3p::LapackLU<cla3p::dns::RdMatrix> luSolver;

	/*
	 * Decompose A into LU product.
	 */
	luSolver.decompose(A);

	{
		/*
		 * Single column (vector) rhs.
		 * Overwrite x with the solution (A^{-1} * b).
		 */

		cla3p::dns::RdVector x = b;

		luSolver.solve(x);

		std::cout << "Dense Vector rhs::Absolute Error: "
			      << (b - A * x).evaluate().normOne() << std::endl;
	}

	{
		/*
		 * Multiple column (matrix) rhs.
		 * Overwrite X with the solution (A^{-1} * B).
		 */

		cla3p::dns::RdMatrix X = B;

		luSolver.solve(X);

		std::cout << "Dense Matrix rhs::Absolute Error: "
			      << (B - A * X).evaluate().normOne() << std::endl;
	}

	return 0;
}
