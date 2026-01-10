/**
 * @example ex01a_dense_vector_create.cpp
 */

#include <iostream>
#include "culite/dense.hpp"

int main()
{
	/*
	 * Double precision real empty vector
	 */

	culite::dns::RdVector X;

	std::cout << X.info("X");

	/*
	 * Single precision real vector with size 3 (uninitialized values)
	 */

	culite::dns::RfVector Y(3);

	std::cout << Y.info("Y");

	/*
	 * Allocate space for X (size 5, uninitialized values)
	 */

	X = culite::dns::RdVector(5);

	std::cout << X.info("X");

	return 0;
}
