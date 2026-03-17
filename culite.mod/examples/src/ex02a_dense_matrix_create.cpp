/**
 * @example ex02a_dense_matrix_create.cpp
 */

#include <iostream>
#include <culite/dense.hpp>

int main()
{
    /*
     * Double precision real empty matrix
     */
    culite::dns::RdMatrix A;
    std::cout << A.info("A");

    /*
     * (3x4) single precision real matrix (uninitialized values)
     */
    culite::dns::RfMatrix B(3,4);
    std::cout << B.info("B");

    /*
     * Allocate space for A (5x2, uninitialized values)
     */
    A = culite::dns::RdMatrix(5,2);
    std::cout << A.info("A");

    return 0;
}
