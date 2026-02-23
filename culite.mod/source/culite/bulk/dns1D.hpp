/*
 * Copyright (c) 2025-2026 Simulisoft
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CULITE_BULK_DNS1D_HPP_
#define CULITE_BULK_DNS1D_HPP_

/**
 * @file
 */

#include "culite/types/integer.hpp"
#include "culite/types/scalar.hpp"
#include "culite/support/utils.hpp"
#include "culite/bulk/dns1D_impl.hpp"
#include "culite/generic/cublas_handler.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/

template <typename T_Scalar>
void copy1D(int_t n, const T_Scalar *x, T_Scalar *b)
{
	if(n > 0) {
		memCopyD2D(n, x, b);
	}
}
/*-------------------------------------------------*/
//
// Norm 1
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normOne1D(int_t n, const T_Scalar *x)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	globalCuBlasHandler().asum<T_Scalar>(n, x, 1, &ret);
	return ret;
}
/*-------------------------------------------------*/
//
// Norm Inf (CAUTION: calls cudaMemcpy to move result from device to host)
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normInf1D(int_t n, const T_Scalar *x)
{
	int_t idx = 0;
	globalCuBlasHandler().iamax<T_Scalar>(n, x, 1, &idx);
	T_Scalar ret = makeScalar<T_Scalar>(0);
	memCopyD2H(1, x + idx - 1, &ret);
	return arith::abs(ret);
}
/*-------------------------------------------------*/
//
// Norm Euclidean
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normEuc1D(int_t n, const T_Scalar *x)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	globalCuBlasHandler().nrm2<T_Scalar>(n, x, 1, &ret);
	return ret;
}
/*-------------------------------------------------*/
//
// Scale 1D
//
template <typename T_Scalar>
void scale1D(int_t n, const T_Scalar& alpha, T_Scalar *x)
{
	globalCuBlasHandler().scal<T_Scalar>(n, &alpha, x, 1);
}
/*-------------------------------------------------*/
//
// Conjugate 1D
//
template <typename T_Scalar>
void conjugate1D(int_t n, T_Scalar* x)
{
	launch_conjugate_kernel_1d(n, x);
}
/*-------------------------------------------------*/
//
// get Real 1D
//
template <typename T_Scalar>
void getReal1D(int_t n, const T_Scalar* x, typename TypeTraits<T_Scalar>::real_type* y)
{
	launch_get_real_kernel_1d<T_Scalar>(n, x, y);
}
/*-------------------------------------------------*/
//
// get Imag 1D
//
template <typename T_Scalar>
void getImag1D(int_t n, const T_Scalar* x, typename TypeTraits<T_Scalar>::real_type* y)
{
	launch_get_imag_kernel_1d<T_Scalar>(n, x, y);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void geevCalculateComplexEigenvalues(int_t n, 
                                     const typename TypeTraits<T_Scalar>::real_type* wri, 
                                     T_Scalar *w)
{
    launch_geev_calculate_complex_eigenvalues_kernel(n, wri, w);
}
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS1D_HPP_
