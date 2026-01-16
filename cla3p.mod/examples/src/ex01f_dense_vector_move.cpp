/**
 * @example ex01f_dense_vector_move.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
	cla3p::dns::RfVector x(3);
	x = 1;

	std::cout << x.info("x") << x;

	/*
	 * Move x to y
	 */

	cla3p::dns::RfVector y = x.move();

	std::cout << x.info("x") << x;
	std::cout << y.info("y") << y;

	/*
	 * Move y to z
	 * z is non-empty with non matching dimension, so operation will fail
	 */

	cla3p::dns::RfVector z(2);

	try {
		z = y.move();
	} catch (...) {
		std::cout << "Invalid move exception caught" << std::endl;
	}

	return 0;
}
