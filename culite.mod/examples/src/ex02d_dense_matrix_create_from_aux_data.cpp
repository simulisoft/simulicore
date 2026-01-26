/**
 * @example ex02d_dense_matrix_create_from_aux_data.cpp
 */

#include <iostream>
#include <cla3p/support.hpp>
#include <culite/support.hpp>
#include <culite/dense.hpp>

int main()
{
	/*
	 * Allocate space for a & b and assume are filled with some values
	 * Values in a and b are assumed to be in column-major order
	 */
	cla3p::uint_t lda = 7;
	cla3p::uint_t ldb = 5;
	cla3p::real_t *a_host = cla3p::i_calloc_t<cla3p::real_t>(lda * 4); 
	cla3p::real_t *b_host = cla3p::i_calloc_t<cla3p::real_t>(ldb * 5); 

	for(cla3p::uint_t j = 0, icnt = 0; j < 4; j++)
		for(cla3p::uint_t i = 0; i < 3; i++)
			a_host[lda * j + i] = icnt++;

	for(cla3p::uint_t j = 0, icnt = 0; j < 5; j++)
		for(cla3p::uint_t i = 0; i < 5; i++)
			b_host[ldb * j + i] = icnt++;

	/*
	 * Copy a & b to device
	 */
	culite::real_t *a_device = culite::device_alloc_t<culite::real_t>(lda * 4);
	culite::real_t *b_device = culite::device_alloc_t<culite::real_t>(ldb * 5);
	culite::memCopyH2D(3, 4, a_host, lda, a_device, lda);
	culite::memCopyH2D(5, 5, b_host, ldb, b_device, ldb);

	/*
	 * Assign pointer a in matrix A but do not bind
	 * A simply hosts a, need to manually dealloc a
	 */
	culite::dns::RdMatrix A = culite::dns::RdMatrix(3, 4, a_device, lda, false);
	std::cout << A.info("A") << A;

	/*
	 * Assign pointer b in matrix B with property and bind
	 * B takes ownership of b, no free call for b is required
	 */
	cla3p::Property prB = cla3p::Property::General();
	culite::dns::RdMatrix B = culite::dns::RdMatrix(5, 5, b_device, ldb, true, prB);
	std::cout << B.info("B") << B;

	/* 
	 * Free a and exit
	 */

	cla3p::i_free(a_host);
	cla3p::i_free(b_host);
	culite::device_free(a_device);
	// b_device is freed by B's destructor

	return 0;
}
