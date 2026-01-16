/**
 * @example ex01j_dense_vector_algebra_add.cpp
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
	 * Perform the operation (x + 2 * y) using operators and the add function respectively
	 */

	cla3p::dns::RdVector z1 = x + 2. * y;
	std::cout << "z1:\n" << z1;

	cla3p::dns::RdVector z2 = cla3p::ops::add(1., x, 2., y);
	std::cout << "z2:\n" << z2 << "\n";

	/*
	 * Perform the operation (z1 += 3 * x) using operators and the update function respectively
	 */

	z1 += 3. * x;
	std::cout << "z1:\n" << z1;

	cla3p::ops::update(3., x, z2);
	std::cout << "z2:\n" << z2;

	return 0;
}
