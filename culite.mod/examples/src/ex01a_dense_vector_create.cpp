/**
 * @example ex01a_dense_vector_create.cpp
 */

#include <iostream>
#include <culite/dense.hpp>

int main()
{
	/*
	 * Double precision real empty vector
	 */

	culite::dns::RdVector x;

	std::cout << x.info("x");

	/*
	 * Single precision real vector with size 3 (uninitialized values)
	 */

	culite::dns::RfVector y(3);

	std::cout << y.info("y");

	/*
	 * Allocate space for x (size 5, uninitialized values)
	 */

	x = culite::dns::RdVector(5);

	std::cout << x.info("x");

	return 0;
}
