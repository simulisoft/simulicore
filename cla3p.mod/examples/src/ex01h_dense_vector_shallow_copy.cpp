/**
 * @example ex01h_dense_vector_shallow_copy.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
	cla3p::dns::RfVector x(3);
	x = 1;

	std::cout << x.info("x") << x;

	/*
	 * Shallow copy x to y (x & y share contents)
	 */

	cla3p::dns::RfVector y = x.rcopy();

	std::cout << y.info("y") << y << "\n";

	y(0) = 0;
	y(2) = 0;

	std::cout << "x:\n" << x;
	std::cout << "y:\n" << y;

	return 0;
}
