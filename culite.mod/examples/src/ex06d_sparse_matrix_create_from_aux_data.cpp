/**
 * @example ex06d_sparse_matrix_create_from_aux_data.cpp
 */

#include <iostream>
#include <culite/support.hpp>
#include <culite/sparse.hpp>

int main()
{
	/*
	 * Allocate space for two csc matrix representations and assume are filled with some values
	 */

	culite::uint_t nr  = 15; // number of rows
	culite::uint_t nc  = 12; // number of columns
	culite::uint_t nnz = 10; // number of non-zeros

	culite::int_t  *rowptrA = culite::device_alloc_t<culite::int_t >(nr + 1); 
	culite::int_t  *colidxA = culite::device_alloc_t<culite::int_t >(nnz); 
	culite::real_t *valuesA = culite::device_alloc_t<culite::real_t>(nnz); 

	culite::int_t  *colptrB = culite::device_alloc_t<culite::int_t >(nc + 1); 
	culite::int_t  *rowidxB = culite::device_alloc_t<culite::int_t >(nnz); 
	culite::real_t *valuesB = culite::device_alloc_t<culite::real_t>(nnz); 

	/*
	 * Assign csr pointers in matrix A but do not bind
	 * A simply hosts, need to manually dealloc csr vectors
	 */

	culite::csr::RdMatrix A = culite::csr::RdMatrix(nr, nc, rowptrA, colidxA, valuesA, false);

	std::cout << A.info("A");

	/*
	 * Assign pointer b in matrix B with property and bind
	 * Assign csc pointers in matrix B with property and bind 
	 * B takes ownership of csc vectors, no free call for the csc vectors is required
	 */

	cla3p::Property prB = cla3p::Property::SymmetricLower();
	culite::csc::RdMatrix B = culite::csc::RdMatrix(nc, nc, colptrB, rowidxB, valuesB, true, prB);

	std::cout << B.info("B");

	/* 
	 * Free a and exit
	 */

	culite::device_free(rowptrA);
	culite::device_free(colidxA);
	culite::device_free(valuesA);

	return 0;
}
