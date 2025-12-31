/**
 * @example ex09a_lowrank_matrix_create.cpp
 */

#include <iostream>
#include "cla3p/dense.hpp"
#include "cla3p/lra.hpp"

int main()
{
	/*
	 * Double precision real empty low-rank matrix
	 */

	cla3p::lra::RdMatrix A;

	std::cout << A.info("A");

	/*
	 * Rank 2 (3x4) single precision real low-rank matrix (uninitialized values)
	 */

	cla3p::lra::XxMatrix<cla3p::dns::RfMatrix> B(3,4,2);

	std::cout << B.info("B");

	/*
	 * Allocate space for A (rank 2, 5x3, uninitialized values)
	 */

	A = cla3p::lra::RdMatrix(5,3,2);

	std::cout << A.info("A");

	return 0;
}
