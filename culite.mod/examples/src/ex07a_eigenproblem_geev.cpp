/**
 * @example ex07a_eigenproblem_geev.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>
#include <culite/eigsol.hpp>

int main()
{
	/*
	 * Create a random (5x5) dense matrix on host and transfer to device
	 */
	const cla3p::dns::RdMatrix HostA = cla3p::dns::RdMatrix::random(5,5);
    const cla3p::dns::RdMatrix HostB = cla3p::dns::RdMatrix::random(7,7);

	culite::dns::RdMatrix A;
    culite::dns::RdMatrix B;
	HostA >> A;
    HostB >> B;
	std::cout << "A:\n" << A;
    std::cout << "B:\n" << B;

	/*
	 * Create an eigensolver that computes right eigenvectors (default)
	 */
	culite::LapackGeev<culite::dns::RdMatrix> geev;

	/*
	 * Decompose A
	 */
	geev.decompose(A);

	std::cout << "Eigenvalues (A):\n" << geev.eigenvalues();
	std::cout << "Right Eigenvectors (A):\n" << geev.rightEigenvectors();

	/*
	 * Create an eigensolver that computes only eigenvalues (no eigenvectors)
	 */
    geev.setCalcLeft(false);
    geev.setCalcRight(false);

	geev.decompose(B);

	std::cout << "Eigenvalues only (B):\n" << geev.eigenvalues();

	return 0;
}
