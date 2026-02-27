/**
 * @example ex01i_dense_vector_algebra_scale.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/algebra.hpp>

int main()
{
	cla3p::dns::RdVector x(3);
	x = 3.;
	std::cout << "x:\n" << x << "\n";

	/*
	 * Scale x using operators and the scale function respectively
	 */
	x *= 2.;
	std::cout << "x *= 2:\n" << x << "\n";

	x.iscale(.5);
	std::cout << "x.iscale(.5):\n" << x << "\n";

	cla3p::dns::RdVector y = 2. * x;
	std::cout << "y:\n" << y;

	return 0;
}
