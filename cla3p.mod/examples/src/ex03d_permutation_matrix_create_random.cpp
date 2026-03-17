/**
 * @example ex03d_permutation_matrix_create_random.cpp
 */

#include <iostream>
#include <cla3p/perms.hpp>
#include <cla3p/support/imalloc.hpp>
#include <cla3p/bulk/dns_io.hpp>

/*-----------------------------------------------------*/
cla3p::int_t* permToDense(const cla3p::prm::PiMatrix& P)
{
    cla3p::int_t n = P.size();
    cla3p::int_t *ret = cla3p::i_calloc_t<cla3p::int_t>(n * n);
    for(cla3p::int_t i = 0; i < n; i++)
        ret[i * n + P(i)] = 1;
    return ret;
}
/*-----------------------------------------------------*/
int main()
{
    /* 
     * (4x4) permutation matrix 
     * with randomly generated unique values
     */

    cla3p::prm::PiMatrix P = cla3p::prm::PiMatrix::random(4);

    std::cout << P.info("P") << P;

    cla3p::int_t *P4x4 = permToDense(P);
    std::cout << "As a dense matrix:\n";
    cla3p::blk::dns::print_to_stream(std::cout, cla3p::uplo_t::Full, 4, 4, P4x4, 4);
    cla3p::i_free(P4x4);

    return 0;
}
