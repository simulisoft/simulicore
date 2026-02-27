/**
 * @example ex06m_sparse_matrix_algebra_vmult.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/sparse.hpp>
#include <culite/dense.hpp>
#include <culite/sparse.hpp>
#include <culite/algebra.hpp>

int main()
{
    cla3p::coo::RdMatrix hostAcoo(5, 5);
    hostAcoo.insert(0,0,1.0);
    hostAcoo.insert(1,1,2.0);
    hostAcoo.insert(2,1,3.0);
    hostAcoo.insert(1,3,4.0);
    hostAcoo.insert(0,0,5.0);
    cla3p::csr::RdMatrix hostA = hostAcoo.toCsr();

	cla3p::dns::RdVector hostX(5);
	hostX = 2.;

	culite::csr::RdMatrix A;
    hostA >> A; // Transfer to GPU
	std::cout << "A:\n" << A << "\n";

	culite::dns::RdVector x(5);
	hostX >> x; // Transfer to GPU
	std::cout << "x:\n" << x << "\n";

	/*
	 * Perform the operation (A * x) using operators and the mult function respectively
	 */
	culite::dns::RdVector y1 = A * x;
	std::cout << "y1:\n" << y1 << "\n";

	culite::dns::RdVector y2(5);
	culite::ops::mult(1., cla3p::op_t::N, A, x, 0., y2);
	std::cout << "y2:\n" << y2 << "\n";

	/*
	 * Perform the operation (y1 += A * x) using operators and the mult function respectively
	 */
	y1 += A * x;
	std::cout << "y1:\n" << y1 << "\n";

	culite::ops::mult(1., cla3p::op_t::N, A, x, 1., y2);
	std::cout << "y2:\n" << y2 << "\n";

	return 0;
}
