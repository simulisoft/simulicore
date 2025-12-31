/**
 * @example ex01i_dense_vector_algebra_scale.cpp
 */

#include <iostream>
#include "cla3p/dense.hpp"
#include "cla3p/algebra.hpp"

int main()
{
	cla3p::dns::RdVector X(3);

	X = 3.;
	std::cout << "X:\n" << X << "\n";

	/*
	 * Scale X using operators and the scale function respectively
	 */

	X *= 2.;
	std::cout << "X *= 2:\n" << X << "\n";

	X.iscale(.5);
	std::cout << "X.iscale(.5):\n" << X << "\n";

	cla3p::dns::RdVector Y = 2. * X;

	std::cout << "Y:\n" << Y;

	return 0;
}
