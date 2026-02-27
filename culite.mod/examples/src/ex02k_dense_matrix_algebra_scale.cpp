/**
 * @example ex02k_dense_matrix_algebra_scale.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>
#include <culite/algebra.hpp>

int main()
{
	cla3p::dns::RdMatrix hostA(3,3);
	hostA = 3.;
    culite::dns::RdMatrix A;
    hostA >> A; // Transfer to GPU
	std::cout << "A:\n" << A << "\n";

	/*
	 * Scale A using operators and the scale function respectively
	 */
	A *= 2.;
	std::cout << "A *= 2:\n" << A;

	A.iscale(.5);
	std::cout << "A.iscale(.5):\n" << A << "\n";

	culite::dns::RdMatrix B = 2. * A ;
	std::cout << "B:\n" << B;

	return 0;
}
