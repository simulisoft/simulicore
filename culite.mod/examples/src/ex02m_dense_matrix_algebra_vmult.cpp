/**
 * @example ex02m_dense_matrix_algebra_vmult.cpp
 */

#include <iostream>
#include <cla3p/dense.hpp>
#include <culite/dense.hpp>
#include <culite/algebra.hpp>

int main()
{
    cla3p::dns::RdMatrix hostA(3, 3);
    cla3p::dns::RdVector hostX(3);
    hostA = 3.;
    hostX = 2.;

    culite::dns::RdMatrix A;
    culite::dns::RdVector x;
    hostA >> A; // Transfer to GPU
    hostX >> x; // Transfer to GPU

    std::cout << "A:\n" << A;
    std::cout << "x:\n" << x << "\n";

    /*
     * Perform the operation (A * x) using operators and the mult function respectively
     */
    culite::dns::RdVector y1 = A * x;
    std::cout << "y1:\n" << y1;

    culite::dns::RdVector y2(3);
    culite::ops::mult(1., cla3p::op_t::N, A, x, 0., y2);
    std::cout << "y2:\n" << y2 << "\n";

    /*
     * Perform the operation (y1 += A * x) using operators and the mult function respectively
     */
    y1 += A * x;
    std::cout << "y1:\n" << y1;

    culite::ops::mult(1., cla3p::op_t::N, A, x, 1., y2);
    std::cout << "y2:\n" << y2;

    return 0;
}
