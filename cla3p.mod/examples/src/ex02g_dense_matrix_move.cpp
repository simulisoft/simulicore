/**
 * @example ex02g_dense_matrix_move.cpp
 */

#include <iostream>
#include "cla3p/dense.hpp"

int main()
{
	cla3p::dns::RfMatrix A(3,3);
	A = 2;

	std::cout << A.info("A") << A;

	/*
	 * Move A to B
	 */

	cla3p::dns::RfMatrix B = A.move();

	std::cout << A.info("A") << A;
	std::cout << B.info("B") << B;

	/*
	 * Move B to C
	 * C is non-empty with non matching dimension, so operation will fail
	 */

	cla3p::dns::RfMatrix C(2,2);

	try {
		C = B.move();
	} catch(...) {
		std::cout << "Invalid move exception caught" << std::endl;
	}

	return 0;
}
