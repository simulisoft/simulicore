/**
 * @example ex01c_dense_vector_create_from_aux_data.cpp 
 */

#include <iostream>
#include <cla3p/support.hpp>
#include <culite/support.hpp>
#include <culite/dense.hpp>

int main()
{
	/*
	 * Allocate space for x & y on host
	 */

	cla3p::real_t *x_host = cla3p::i_calloc_t<cla3p::real_t>(5); 
	cla3p::real_t *y_host = cla3p::i_calloc_t<cla3p::real_t>(3); 

	for(cla3p::uint_t i = 0; i < 5; i++)
		x_host[i] = i;

	for(cla3p::uint_t i = 0; i < 3; i++)
		y_host[i] = 3 - i;

	/*
	 * Allocate space for x & y on device and copy data from host to device
	 */
	culite::real_t *x_device = culite::device_alloc_t<culite::real_t>(5);
	culite::real_t *y_device = culite::device_alloc_t<culite::real_t>(3);
	culite::memCopyH2D(5, x_host, x_device);
	culite::memCopyH2D(3, y_host, y_device);

	/*
	 * Assign pointer x in vector x but do not bind
	 * x simply hosts x, need to manually dealloc x
	 */

	culite::dns::RdVector x(5, x_device, false);
	std::cout << x.info("x") << x;

	/*
	 * Assign pointer y in vector y and bind
	 * y takes ownership of y, no free call for y is required
	 */

	culite::dns::RdVector y(3, y_device, true);
	std::cout << y.info("y") << y;

	/* 
	 * Free x and exit
	 */

	cla3p::i_free(x_host);
	cla3p::i_free(y_host);
	culite::device_free(x_device);
	// y_device is freed by y's destructor
    
	return 0;
}
