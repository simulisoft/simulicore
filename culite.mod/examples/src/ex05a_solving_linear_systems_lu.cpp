/**
 * @example ex05a_solving_linear_systems_lu.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/algebra.hpp>
#include <culite/dense.hpp>
#include <culite/linsol.hpp>

int main()
{
	/*
	 * Create a random dense objects on host
	 */
	const cla3p::dns::RdMatrix HostA  = cla3p::dns::RdMatrix::random(5,5);
	const cla3p::dns::RdVector HostB1 = cla3p::dns::RdVector::random(5);
	const cla3p::dns::RdMatrix HostB2 = cla3p::dns::RdMatrix::random(5,3);

	/*
	 * Transfer to device
	 */
	culite::dns::RdMatrix A ;
	culite::dns::RdVector B1;
	culite::dns::RdMatrix B2;
	HostA >> A;
	HostB1 >> B1;
	HostB2 >> B2;

	/*
	 * Instantiate LU solver
	 */
	culite::LapackLU<culite::dns::RdMatrix> luSolver;

	/*
	 * Decompose A into LU product
	 */
	luSolver.decompose(A);

	{
		/*
		 * Single column (vector) rhs
		 * Overwrite X with the solution (A^{-1} * B1)
		 */
		culite::dns::RdVector X = B1;
		luSolver.solve(X);
		// Check error norm on host
		cla3p::dns::RdVector HostX;
		X >> HostX;
		std::cout << "Dense Vector rhs::Absolute Error: "
			<< (HostB1 - HostA * HostX).evaluate().normOne() << std::endl;
	}

	{
		/*
		 * Multiple column (matrix) rhs
		 * Overwrite X with the solution (A^{-1} * B2)
		 */
		culite::dns::RdMatrix X = B2;
		luSolver.solve(X);
		// Check error norm on host
		cla3p::dns::RdMatrix HostX;
		X >> HostX;
		std::cout << "Dense Matrix rhs::Absolute Error: "
			<< (HostB2 - HostA * HostX).evaluate().normOne() << std::endl;
	}

	return 0;
}
