/**
 * @example ex02l_dense_matrix_algebra_add.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>
#include <culite/algebra.hpp>

int main()
{
	cla3p::dns::RdMatrix hostA(3, 3);
	cla3p::dns::RdMatrix hostB(3, 3);
	culite::dns::RdMatrix A;
	culite::dns::RdMatrix B;

	hostA = 3.;
	hostB = 2.;
	hostA >> A; // Transfer to GPU
	hostB >> B; // Transfer to GPU
	std::cout << "A:\n" << A;
	std::cout << "B:\n" << B << "\n";

	/*
	 * Perform the operation (A + 2 * B) using operators and the add function respectively
	 */
	culite::dns::RdMatrix C1 = A + 2. * B;
	std::cout << "C1:\n" << C1;

	culite::dns::RdMatrix C2(3, 3);
	culite::ops::add(cla3p::op_t::N, 1., A, 
                     cla3p::op_t::N, 2., B, C2);
	std::cout << "C2:\n" << C2 << "\n";

	/*
	 * Perform the operation (Cx += 3 * A) using operators and the update function respectively
	 */
	C1 += 3. * A;
	std::cout << "C1:\n" << C1;

	culite::ops::update(3., A, C2);
	std::cout << "C2:\n" << C2;

	return 0;
}
