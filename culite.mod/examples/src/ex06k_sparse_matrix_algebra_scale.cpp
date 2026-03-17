/**
 * @example ex06k_sparse_matrix_algebra_scale.cpp
 */

#include <iostream>
#include <cla3p/sparse.hpp>
#include <culite/sparse.hpp>
#include <culite/algebra.hpp>

int main()
{
    cla3p::coo::RdMatrix hostAcoo(5, 5);
    hostAcoo.insert(0,0,1.0);
    hostAcoo.insert(1,1,2.0);
    hostAcoo.insert(2,1,3.0);
    hostAcoo.insert(1,3,4.0);
    hostAcoo.insert(0,0,5.0);
    cla3p::csr::RdMatrix hostA = hostAcoo.toCsr();

    culite::csr::RdMatrix A;
    hostA >> A; // Transfer to GPU
    std::cout << "A:\n" << A << "\n";

    /*
     * Scale A using operators and the scale function respectively
     */
    A *= 2.;
    std::cout << "A *= 2:\n" << A << "\n";

    A.iscale(.5);
    std::cout << "A.iscale(.5):\n" << A << "\n";

    culite::csr::RdMatrix B = 2. * A ;
    std::cout << "B:\n" << B << "\n";

    return 0;
}
