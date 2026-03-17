/**
 * @example ex08a_sv_decomposition.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>
#include <culite/svd.hpp>

int main()
{
	/*
	 * Create a random (5x3) dense matrix on host and transfer to device
	 */
	const cla3p::dns::RdMatrix HostA = cla3p::dns::RdMatrix::random(5,3);

	culite::dns::RdMatrix A;
	HostA >> A;
	std::cout << "A:\n" << A;

	/*
	 * Set a simple svd object for A
	 */
	culite::DefaultSVD<culite::dns::RdMatrix> svdA;

	/*
	 * Decompose A into USV* product
	 */
	svdA.decompose(A);

	std::cout << "SVD Decomposition of A\n";
	std::cout << "----------------------\n";
	std::cout
		<< "S:\n" << svdA.singularValues()
		<< "U:\n" << svdA.leftSingularVectors()
		<< "V:\n" << svdA.rightSingularVectors();

	/*
	 * Set a svd object with no calculation for V
	 */
	culite::svdPolicy_t policyU = culite::svdPolicy_t::Limited;
	culite::svdPolicy_t policyV = culite::svdPolicy_t::NoCalculation;
	culite::DefaultSVD<culite::dns::RdMatrix> svdB(policyU, policyV);

	svdB.decompose(A);

	std::cout << "SVD Decomposition of A (no V)\n";
	std::cout << "-----------------------------\n";
	std::cout
		<< "S:\n" << svdB.singularValues()
		<< "U:\n" << svdB.leftSingularVectors()
		<< "V(empty):\n" << svdB.rightSingularVectors();

	return 0;
}
