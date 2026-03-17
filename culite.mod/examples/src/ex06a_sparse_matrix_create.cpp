/**
 * @example ex06a_sparse_matrix_create.cpp
 */

#include <iostream>
#include <culite/sparse.hpp>

int main()
{
    /*
     * Double precision real empty matrix
     */

    culite::csr::RdMatrix A;
    culite::csc::RdMatrix B;

    std::cout << A.info("A");
    std::cout << B.info("B");

    /*
     * (3x4) single precision real matrix
     */

    culite::csr::RfMatrix C(3,4,7);

    std::cout << C.info("C");

    /*
     * Initialize a (5x2) empty sparse matrix with space for 7 non-zeros
     */

    B = culite::csc::RdMatrix(5,2,7);

    std::cout << B.info("B");

    return 0;
}
