/**
 * @example ex01c_dense_vector_create_from_aux_data.cpp 
 */

#include <iostream>
#include <cla3p/support.hpp>
#include <cla3p/dense.hpp>

int main()
{
	/*
	 * Allocate space for x & y
	 */
	cla3p::real_t *p_x = cla3p::i_calloc_t<cla3p::real_t>(5); 
	cla3p::real_t *p_y = cla3p::i_calloc_t<cla3p::real_t>(3); 

	for(cla3p::uint_t i = 0; i < 5; i++)
		p_x[i] = i;
	for(cla3p::uint_t i = 0; i < 3; i++)
		p_y[i] = 3 - i;

	/*
	 * Assign pointer p_x to vector x but do not bind.
	 * (bind = false) x does not take ownership of p_x, need to manually deallocate p_x.
	 */
	cla3p::dns::RdVector x(5, p_x, false);
	std::cout << x.info("x") << x;

	/*
	 * Assign pointer p_y to vector y with bind.
	 * (bind = true) y takes ownership of p_y, no free call for p_y is required.
	 */
	cla3p::dns::RdVector y(3, p_y, true);
	std::cout << y.info("y") << y;

	/*
	 * Free p_x and exit.
	 */
	cla3p::i_free(p_x);

	return 0;
}
