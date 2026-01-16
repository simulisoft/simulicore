/**
 * @example ex01g_dense_vector_block_operations.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>

int main()
{
	cla3p::dns::RdVector x(5);

	for(cla3p::uint_t i = 0; i < 5; i++)
		x(i) = i;

	std::cout << "x:\n" << x << "\n";

	/*
	 * Get a copy of x[1:3) in xb
	 * Get a reference of x[2:5) in xr
	 */
	cla3p::dns::RdVector xb = x.block(1, 2);
	cla3p::dns::RdVector xr = x.rblock(2, 3);
	std::cout << "xb:\n" << xb;
	std::cout << "xr:\n" << xr;
	/*
	 * Get a guarded reference of xref[2:5) in xg
	 */
	const cla3p::dns::RdVector& xref = x;
	cla3p::Guard<cla3p::dns::RdVector> xg = xref.rblock(2, 3);
	std::cout << "xg:\n" << xg.get() << "\n";

	/*
	 * Change values in blocks
	 */
	xr = -1;
	std::cout << "x:\n" << x;

	/*
	 * Set values of xb to x, starting at index 1
	 */
	xb = -2;
	x.setBlock(1, xb);
	std::cout << "x:\n" << x;

	/*
	 * Set values of xb to x, starting at index 2
	 */
	xb = -3;
	x.rblock(2, 2) = xb;
	std::cout << "x:\n" << x;

	return 0;
}
