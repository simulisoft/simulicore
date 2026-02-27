/**
 * @example ex01k_dense_vector_algebra_dot.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/algebra.hpp>

int main()
{
	cla3p::dns::RdVector x(3);
	cla3p::dns::RdVector y(3);

	x = 3.;
	y = 2.;

	std::cout << "x:\n" << x;
	std::cout << "y:\n" << y << "\n";

	/*
	 * Perform the operation (x' * y) using operators and the dot function respectively
	 */
	cla3p::real_t d1 = x.transpose() * y;
	std::cout << "d1: " << d1 << std::endl;

	cla3p::real_t d2 = cla3p::ops::dot(x, y);
	std::cout << "d2: " << d2 << std::endl;

	return 0;
}
