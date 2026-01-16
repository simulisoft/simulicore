/**
 * @example ex04a_guard_const_vector.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
	/*
	 * Create a constant vector
	 */
	const cla3p::dns::RdVector x = cla3p::dns::RdVector::random(5);

	/*
	 * Get a shallow copy of the constant vector using the Guard class
	 */
	cla3p::Guard<cla3p::dns::RdVector> y_grd = x.rcopy();

	/*
	 * Get an immutable reference of the copied vector
	 */
	const cla3p::dns::RdVector& y = y_grd.get();

	/*
	 * Get info of the two vectors
	 */
	std::cout << x.info("x") << y.info("y") << std::endl;

	return 0;
}
