/**
 * @example ex06d_sparse_matrix_create_from_aux_data.cpp
 */

#include <iostream>
#include <cla3p/support.hpp>
#include <cla3p/sparse.hpp>

int main()
{
    /*
     * Allocate space for two CSC matrix representations and assume they are filled with some values.
     */
    cla3p::uint_t nr  = 15; // number of rows
    cla3p::uint_t nc  = 12; // number of columns
    cla3p::uint_t nnz = 10; // number of non-zeros

    cla3p::int_t  *colptrA = cla3p::i_calloc_t<cla3p::int_t >(nc + 1); 
    cla3p::int_t  *rowidxA = cla3p::i_calloc_t<cla3p::int_t >(nnz); 
    cla3p::real_t *valuesA = cla3p::i_calloc_t<cla3p::real_t>(nnz); 

    cla3p::int_t  *colptrB = cla3p::i_calloc_t<cla3p::int_t >(nc + 1); 
    cla3p::int_t  *rowidxB = cla3p::i_calloc_t<cla3p::int_t >(nnz); 
    cla3p::real_t *valuesB = cla3p::i_calloc_t<cla3p::real_t>(nnz); 

    /*
     * Assign CSC pointers to matrix A but do not bind.
     * (bind = false) A does not take ownership, need to manually deallocate CSC vectors.
     */
    cla3p::csc::RdMatrix A(nr, nc, colptrA, rowidxA, valuesA, false);
    std::cout << A.info("A");

    /*
     * Assign CSC pointers to matrix B with property and bind.
     * (bind = true) B takes ownership of CSC vectors, no free call for the CSC vectors is required.
     */
    cla3p::Property prB = cla3p::Property::SymmetricLower();
    cla3p::csc::RdMatrix B(nc, nc, colptrB, rowidxB, valuesB, true, prB);
    std::cout << B.info("B");

    /*
     * Free CSC vectors for A and exit.
     */
    cla3p::i_free(colptrA);
    cla3p::i_free(rowidxA);
    cla3p::i_free(valuesA);

    return 0;
}
