/**
 * @example ex06l_sparse_matrix_algebra_add.cpp
 */

#include <iostream>
#include <cla3p/sparse.hpp>
#include <culite/sparse.hpp>
#include <culite/algebra.hpp>

int main()
{
    cla3p::coo::RdMatrix hostAcoo(5, 5);
    cla3p::coo::RdMatrix hostBcoo(5, 5);

    hostAcoo.insert(0,0,1.0);
    hostAcoo.insert(1,1,2.0);
    hostAcoo.insert(2,1,3.0);
    hostAcoo.insert(1,3,4.0);
    hostAcoo.insert(0,0,5.0);

    hostBcoo.insert(0,0,-1.3);
    hostBcoo.insert(1,2, 2.2);
    hostBcoo.insert(2,1,-3.1);
    hostBcoo.insert(2,4, 4.6);
    hostBcoo.insert(1,1,-5.7);

    cla3p::csr::RdMatrix hostA = hostAcoo.toCsr();
    cla3p::csr::RdMatrix hostB = hostBcoo.toCsr();

	culite::csr::RdMatrix A;
	culite::csr::RdMatrix B;
	hostA >> A; // Transfer to GPU
	hostB >> B; // Transfer to GPU
	std::cout << A.info("A") << A << "\n";
	std::cout << B.info("B") << B << "\n";

	/*
	 * Perform the operation (A + 2 * B) using operators and the add function respectively
	 */
	culite::csr::RdMatrix C1 = 3. * A + 2. * B;
	std::cout << "C1:\n" << C1 << "\n";

	culite::csr::RdMatrix C2 = culite::ops::add(3., A, 2., B);
	std::cout << "C2:\n" << C2 << "\n";

	/*
	 * Perform the operation (Cx += 3 * A) using operators and the update function respectively
	 */
	C1 += 3. * A;
	std::cout << "C1:\n" << C1 << "\n";

	culite::ops::update(3., A, C2);
	std::cout << "C2:\n" << C2 << "\n";

	return 0;
}
