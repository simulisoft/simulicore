/**
 * @example ex02o_dense_matrix_algebra_mtmult.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>
#include <culite/algebra.hpp>

int main()
{
    cla3p::dns::RdMatrix hostA(3, 3);
    cla3p::dns::RdMatrix hostB(3, 3);
    hostA = 3.;
    hostB = 2.;

    culite::dns::RdMatrix A;
    culite::dns::RdMatrix B;
    hostA >> A; // Transfer to GPU
    hostB >> B; // Transfer to GPU
    std::cout << "A:\n" << A;
    std::cout << "B:\n" << B << "\n";

    /*
     * Perform the operation (A' * B) using operators and the mult function respectively
     */
    culite::dns::RdMatrix C1 = A.transpose() * B;
    std::cout << "C1:\n" << C1;

    culite::dns::RdMatrix C2(3, 3);
    culite::ops::mult(1., cla3p::op_t::T, A, cla3p::op_t::N, B, 0., C2);
    std::cout << "C2:\n" << C2;

    return 0;
}
