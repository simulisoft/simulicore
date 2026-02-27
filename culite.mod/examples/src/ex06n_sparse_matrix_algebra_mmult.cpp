/**
 * @example ex06n_sparse_matrix_algebra_mmult.cpp
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
    cla3p::coo::RdMatrix hostBcoo(5, 3);

    hostAcoo.insert(0,0,1.0);
    hostAcoo.insert(1,1,2.0);
    hostAcoo.insert(2,1,3.0);
    hostAcoo.insert(1,3,4.0);
    hostAcoo.insert(0,0,5.0);

    hostBcoo.insert(1,2, 1.3);
    hostBcoo.insert(2,1,-2.1);
    hostBcoo.insert(2,0, 3.5);
    hostBcoo.insert(1,1,-4.4);
    hostBcoo.insert(0,0, 5.3);

    cla3p::csr::RdMatrix hostAcsr = hostAcoo.toCsr();
    cla3p::csr::RdMatrix hostBcsr = hostBcoo.toCsr();

    cla3p::dns::RdMatrix hostBdns(5, 3);
	hostBdns = 2.;

	culite::csr::RdMatrix Acsr;
	culite::csr::RdMatrix Bcsr;
    hostAcsr >> Acsr; // Transfer to GPU
    hostBcsr >> Bcsr; // Transfer to GPU
	std::cout << "A:\n" << Acsr << "\n";
	std::cout << "B:\n" << Bcsr << "\n";

	culite::dns::RdMatrix Bdns(5, 3);
	hostBdns >> Bdns; // Transfer to GPU
	std::cout << "B:\n" << Bdns << "\n";

	/*
	 * Perform the operation (A * B) using operators and the mult function respectively
	 */

	cla3p::op_t opA = cla3p::op_t::N;
	// cla3p::op_t opB = cla3p::op_t::N;

	culite::dns::RdMatrix C1 = Acsr * Bdns;
	std::cout << "C1:\n" << C1 << "\n";

	culite::dns::RdMatrix C2(5, 3);
	culite::ops::mult(1., opA, Acsr, Bdns, 0., C2);
	std::cout << "C2:\n" << C2 << "\n";

    // Sparse-sparse multiplication is currently not supported...
	// culite::csr::RdMatrix C3 = Acsr * Bcsr;
	// std::cout << "C3:\n" << C3 << "\n";
    //
	// culite::csr::RdMatrix C4 = culite::ops::mult(1., opA, Acsr, opB, Bcsr);
	// std::cout << "C4:\n" << C4 << "\n";

	/*
	 * Perform the operation (Cx += A * B) using operators and the mult function respectively
	 */
	C1 += Acsr * Bdns;
	std::cout << "C1:\n" << C1 << "\n";

	culite::ops::mult(1., opA, Acsr, Bdns, 1., C2);
	std::cout << "C2:\n" << C2 << "\n";

    // Sparse-sparse addition is currently not supported...
	// C3 += Acsr * Bcsr;
	// std::cout << "C3:\n" << C3 << "\n";

	return 0;
}
