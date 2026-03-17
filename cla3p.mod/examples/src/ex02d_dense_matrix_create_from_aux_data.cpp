/**
 * @example ex02d_dense_matrix_create_from_aux_data.cpp
 */

#include <iostream>
#include <cla3p/support.hpp>
#include <cla3p/dense.hpp>

int main()
{
    /*
     * Allocate space for a & b and assume they are filled with some values.
     * Values in a and b are assumed to be in column-major order.
     */
    cla3p::uint_t lda = 7;
    cla3p::uint_t ldb = 5;
    cla3p::real_t *a = cla3p::i_calloc_t<cla3p::real_t>(lda * 4); 
    cla3p::real_t *b = cla3p::i_calloc_t<cla3p::real_t>(ldb * 5); 

    for(cla3p::uint_t j = 0, icnt = 0; j < 4; j++)
        for(cla3p::uint_t i = 0; i < 3; i++)
            a[lda * j + i] = icnt++;

    for(cla3p::uint_t j = 0, icnt = 0; j < 5; j++)
        for(cla3p::uint_t i = j; i < 5; i++)
            b[ldb * j + i] = icnt++;

    /*
     * Assign pointer a to matrix A but do not bind.
     * (bind = false) A does not take ownership of a, need to manually deallocate a.
     */
    cla3p::dns::RdMatrix A(3, 4, a, lda, false);
    std::cout << A.info("A") << A;

    /*
     * Assign pointer b to matrix B with property and bind.
     * (bind = true) B takes ownership of b, no free call for b is required.
     */
    cla3p::Property prB = cla3p::Property::SymmetricLower();
    cla3p::dns::RdMatrix B(5, 5, b, ldb, true, prB);
    std::cout << B.info("B") << B;

    /*
     * Free a and exit.
     */
    cla3p::i_free(a);

    return 0;
}
