/**
 * @example ex06m_sparse_matrix_algebra_vmult.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/sparse.hpp>
#include <cla3p/algebra.hpp>

int main()
{
	cla3p::coo::RdMatrix Acoo(5, 5);

	Acoo.insert(0,0,1.0);
	Acoo.insert(1,1,2.0);
	Acoo.insert(2,1,3.0);
	Acoo.insert(1,3,4.0);
	Acoo.insert(0,0,5.0);

	cla3p::csc::RdMatrix A = Acoo.toCsc();
	std::cout << "A:\n" << A << "\n";

	cla3p::dns::RdVector x(5);
	x = 2.;
	std::cout << "x:\n" << x << "\n";

	/*
	 * Perform the operation (A * x) using operators and the mult function respectively
	 */
	cla3p::dns::RdVector y1 = A * x;
	std::cout << "y1:\n" << y1;

	cla3p::dns::RdVector y2(5);
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
