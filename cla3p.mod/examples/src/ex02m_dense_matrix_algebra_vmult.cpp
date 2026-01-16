/**
 * @example ex02m_dense_matrix_algebra_vmult.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/algebra.hpp>

int main()
{
	cla3p::dns::RdMatrix A(3,3);
	cla3p::dns::RdVector x(3);

	A = 3.;
	x = 2.;
	std::cout << "A:\n" << A;
	std::cout << "x:\n" << x << "\n";

	/*
	 * Perform the operation (A * x) using operators and the mult function respectively
	 */

	cla3p::dns::RdVector y1 = A * x;
	std::cout << "y1:\n" << y1;

	cla3p::dns::RdVector y2(3);
	cla3p::ops::mult(1., cla3p::op_t::N, A, x, 0., y2);
	std::cout << "y2:\n" << y2 << "\n";

	/*
	 * Perform the operation (y1 += A * x) using operators and the mult function respectively
	 */

	y1 += A * x;
	std::cout << "y1:\n" << y1;

	cla3p::ops::mult(1., cla3p::op_t::N, A, x, 1., y2);
	std::cout << "y2:\n" << y2;

	return 0;
}
