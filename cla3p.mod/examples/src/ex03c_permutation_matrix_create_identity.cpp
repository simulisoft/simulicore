/**
 * @example ex03c_permutation_matrix_create_identity.cpp
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
	 * (5x5) permutation matrix (identity -> P[i] = i)
	 */

	cla3p::prm::PiMatrix P = cla3p::prm::PiMatrix::identity(5);

	std::cout << P.info("P") << P;

    cla3p::int_t *Q5x5 = permToDense(P);
	std::cout << "As a dense matrix:\n";
    cla3p::blk::dns::print_to_stream(std::cout, cla3p::uplo_t::Full, 5, 5, Q5x5, 5);
    cla3p::i_free(Q5x5);

	return 0;
}
