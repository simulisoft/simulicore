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

#ifndef CULITE_BULK_DNS2D_HPP_
#define CULITE_BULK_DNS2D_HPP_

/**
 * @file
 */

#include <cla3p/types/enums.hpp>

#include "culite/types/integer.hpp"
#include "culite/support/utils.hpp"
#include "culite/error/exceptions.hpp"
#include "culite/bulk/dns2D_impl.hpp"
#include "culite/generic/cublas_handler.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/

//
// Set pointer (column-major)
//
template <typename T_Scalar>
inline T_Scalar* ptrmv(int_t lda, T_Scalar *a, int_t i, int_t j)
{
	return (a + i + j * lda);
}
/*-------------------------------------------------*/
//
// Get value (column-major)
//
template <typename T_Scalar>
inline T_Scalar& entry(int_t lda, T_Scalar *a, int_t i, int_t j)
{
	return *ptrmv(lda,a,i,j);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void copy2D(::cla3p::uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar *b, int_t ldb)
{
	if(m > 0 && n > 0) {
		if(uplo == ::cla3p::uplo_t::Full) {
			memCopyD2D(m, n, a, lda, b, ldb);
		} else {
			err::CudaException("blk::dns::copy2D for Upper/Lower is not implemented yet.");
		}
	}
}
/*-------------------------------------------------*/
//
// Scale
//
template <typename T_Scalar>
void scale2D(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar coeff)
{
	if(uplo == ::cla3p::uplo_t::Full) {
		launch_scale_matrix_kernel<T_Scalar>(m, n, a, lda, coeff);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::scale2D is only implemented for full matrices.");
	}	
}
/*-------------------------------------------------*/
//
// Norm 1
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normOne2D(::cla3p::prop_t ptype, 
                                                   ::cla3p::uplo_t uplo,
                                                   int_t m, int_t n, 
                                                   const T_Scalar* a, 
                                                   int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_1_norm(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normOne2D is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
//
// Norm Inf
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normInf2D(::cla3p::prop_t ptype, 
												  ::cla3p::uplo_t uplo,
												  int_t m, int_t n, 
												  const T_Scalar* a, 
												  int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_inf_norm(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normInf2D is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
//
// Norm Max
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normMax2D(::cla3p::prop_t ptype, 
	                                               ::cla3p::uplo_t uplo,
                                                   int_t m, int_t n, 
                                                   const T_Scalar* a, 
                                                   int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_max_norm(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normMax2D is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
//
// Norm Frobenius
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normFro2D(::cla3p::prop_t ptype, 
                                                   ::cla3p::uplo_t uplo,
                                                   int_t m, int_t n, 
                                                   const T_Scalar* a, 
                                                   int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_fro_norm(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normFro2D is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
//
// Conjugate 2D
//
template <typename T_Scalar>
//void conjugate2D(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda)
void conjugate2D(::cla3p::uplo_t, int_t, int_t, T_Scalar*, int_t)
{
	// TODO: implement
	throw err::CudaException("blk::dns::conjugate2D is not implemented yet.");
}
/*-------------------------------------------------*/
//
// Transpose 2D
//
template <typename T_Scalar>
void transpose2D(int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar* b, int_t ldb)
{
	T_Scalar alpha = makeScalar<T_Scalar>(1);
	T_Scalar beta  = makeScalar<T_Scalar>(0);
	globalCuBlasHandler().geam(::cla3p::op_t::T,
					           ::cla3p::op_t::N,
								m, n,
					            &alpha, a, lda,
					            &beta,  static_cast<const T_Scalar*>(nullptr), 1,
					            b, ldb);
}
/*-------------------------------------------------*/
//
// Transpose 2D
//
template <typename T_Scalar>
void ctranspose2D(int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar* b, int_t ldb)
{
	T_Scalar alpha = makeScalar<T_Scalar>(1);
	T_Scalar beta  = makeScalar<T_Scalar>(0);
	globalCuBlasHandler().geam<T_Scalar>(::cla3p::op_t::C,
					                     ::cla3p::op_t::N,
								         m, n,
					                     &alpha, a, lda,
					                     &beta,  static_cast<const T_Scalar*>(nullptr), 1,
					                     b, ldb);
}
/*-------------------------------------------------*/
//
// get Real 2D
//
template <typename T_Scalar>
void getReal2D(::cla3p::uplo_t, int_t, int_t, const T_Scalar*, int_t, typename TypeTraits<T_Scalar>::real_type*, int_t)
{
	// TODO: implement
	throw err::CudaException("blk::dns::getReal2D is not implemented yet.");
}
/*-------------------------------------------------*/
//
// get Imag 2D
//
template <typename T_Scalar>
void getImag2D(::cla3p::uplo_t, int_t, int_t, const T_Scalar*, int_t, typename TypeTraits<T_Scalar>::real_type*, int_t)
{
	// TODO: implement
	throw err::CudaException("blk::dns::getImag2D is not implemented yet.");
}
/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS2D_HPP_
