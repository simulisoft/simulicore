/**
 * @example ex10a_geev_decomposition.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <cla3p/eigsol.hpp>

int main()
{
    /*
     * Create a random 4x4 general (non-symmetric) matrix.
     */
    const cla3p::dns::RdMatrix A = cla3p::dns::RdMatrix::random(4, 4);

    std::cout << "A:\n" << A;

    /*
     * Set a simple geev object (right eigenvectors only, the default).
     */
    cla3p::LapackGeev<cla3p::dns::RdMatrix> geev;

    /*
     * Decompose A: compute eigenvalues and right eigenvectors.
     */
    geev.decompose(A);

    std::cout << "GEEV Decomposition of A (right eigenvectors only)\n";
    std::cout << "-------------------------------------------------\n";
    std::cout << "Eigenvalues:\n"      << geev.eigenvalues()      << "\n";
    std::cout << "Right Eigenvectors:\n" << geev.rightEigenvectors() << "\n";

    /*
     * Set a geev object requesting both left and right eigenvectors,
     * with buffer pre-allocation.
     */
    cla3p::LapackGeev<cla3p::dns::RdMatrix> geevLR(true, true);
    geevLR.reserve(4);

    geevLR.decompose(A);

    std::cout << "GEEV Decomposition of A (left and right eigenvectors)\n";
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Eigenvalues:\n"       << geevLR.eigenvalues()       << "\n";
    std::cout << "Left Eigenvectors:\n"  << geevLR.leftEigenvectors()  << "\n";
    std::cout << "Right Eigenvectors:\n" << geevLR.rightEigenvectors() << "\n";

    /*
     * Set a geev object computing eigenvalues only (no eigenvectors).
     */
    cla3p::LapackGeev<cla3p::dns::RdMatrix> geevEigsOnly(false, false);

    geevEigsOnly.decompose(A);

    std::cout << "GEEV Decomposition of A (eigenvalues only)\n";
    std::cout << "------------------------------------------\n";
    std::cout << "Eigenvalues:\n" << geevEigsOnly.eigenvalues() << "\n";

    return 0;
}
