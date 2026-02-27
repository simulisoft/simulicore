/**
 * @example ex01a_dense_vector_create.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
	/*
	 * Double precision real empty vector
	 */
	cla3p::dns::RdVector x;
	std::cout << x.info("x");

	/*
	 * Single precision real vector with size 3 (uninitialized values)
	 */
	cla3p::dns::RfVector y(3);
	std::cout << y.info("y");

	/*
	 * Allocate space for x (size 5, uninitialized values)
	 */
	x = cla3p::dns::RdVector(5);
	std::cout << x.info("x");

	return 0;
}
