/**
 * @example ex01f_dense_vector_move.cpp
 */

#include <iostream>
#include "cla3p/dense.hpp"

int main()
{
	cla3p::dns::RfVector X(3);
	X = 1;

	std::cout << X.info("X") << X;

	/*
	 * Move X to Y
	 */

	cla3p::dns::RfVector Y = X.move();

	std::cout << X.info("X") << X;
	std::cout << Y.info("Y") << Y;

	/*
	 * Move Y to Z
	 * Z is non-empty with non matching dimension, so operation will fail
	 */

	cla3p::dns::RfVector Z(2);

	try {
		Z = Y.move();
	} catch (...) {
		std::cout << "Invalid move exception caught" << std::endl;
	}

	return 0;
}
