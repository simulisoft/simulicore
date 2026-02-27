/**
 * @example ex01e_dense_vector_copy.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
	cla3p::dns::RfVector x(3);
	x = 1;

	std::cout << x.info("x") << x;

	/*
	 * Copy x to y using assignment.
	 * Copy x to z using copy().
	 */
	cla3p::dns::RfVector y = x;
	cla3p::dns::RfVector z = x.copy();
	std::cout << y.info("y") << y;
	std::cout << z.info("z") << z;

	/*
	 * Copy y to z using assignment.
	 * z is already allocated.
	 * No re-allocation performed.
	 * Needs matching sizes.
	 */
	z = y;
	std::cout << z.info("z") << z;

	return 0;
}
