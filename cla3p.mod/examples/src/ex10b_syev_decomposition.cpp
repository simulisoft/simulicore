/**
 * @example ex10b_syev_decomposition.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/eigsol.hpp>
#include <cla3p/algebra.hpp>

int main()
{
    /*
     * Create a random 4x4 symmetric positive-definite matrix: A = B' * B.
     */
    const cla3p::dns::RdMatrix A = cla3p::dns::RdMatrix::random(4, 4, cla3p::Property::SymmetricLower());
    std::cout << A.info("A") << "\n" << A;

    /*
     * Set a simple syev object (eigenvectors are computed by default).
     */
    cla3p::LapackSyev<cla3p::dns::RdMatrix> syev;

    /*
     * Decompose A: compute eigenvalues and eigenvectors.
     */
    syev.decompose(A);

    std::cout << "SYEV Decomposition of A (eigenvalues + eigenvectors)\n";
    std::cout << "----------------------------------------------------\n";
    std::cout << "Eigenvalues (real, ascending):\n" << syev.eigenvalues()  << "\n";
    std::cout << "Eigenvectors:\n"                  << syev.eigenvectors() << "\n";

    /*
     * Set a syev object with buffer pre-allocation and eigenvalues only.
     */
    cla3p::LapackSyev<cla3p::dns::RdMatrix> syevEigsOnly(false);
    syevEigsOnly.reserve(4);
    syevEigsOnly.decompose(A);

    std::cout << "SYEV Decomposition of A (eigenvalues only)\n";
    std::cout << "------------------------------------------\n";
    std::cout << "Eigenvalues (real, ascending):\n" << syevEigsOnly.eigenvalues() << "\n";

    return 0;
}
