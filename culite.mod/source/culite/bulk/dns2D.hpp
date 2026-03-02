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
#include "culite/types/scalar.hpp"
#include "culite/support/utils.hpp"
#include "culite/error/exceptions.hpp"
#include "culite/bulk/dns2D_impl.hpp"
#include "culite/generic/cublas_handler.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace blk {
namespace dns {
/*-------------------------------------------------*/
template <typename T_Scalar>
inline T_Scalar* ptrmv(int_t lda, T_Scalar *a, int_t i, int_t j)
{
	return (a + i + j * lda);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
inline const T_Scalar* ptrmv(int_t lda, const T_Scalar *a, int_t i, int_t j)
{
	return (a + i + j * lda);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
inline T_Scalar& entry(int_t lda, T_Scalar *a, int_t i, int_t j)
{
	return *ptrmv(lda, a, i, j);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
inline const T_Scalar& entry(int_t lda, const T_Scalar *a, int_t i, int_t j)
{
	return *ptrmv(lda, a, i, j);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void copy(::cla3p::uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar *b, int_t ldb)
{
	if(m > 0 && n > 0) {
		if(uplo == ::cla3p::uplo_t::Full) {
			memCopyD2D(m, n, a, lda, b, ldb);
		} else {
            launch_copy_kernel(uplo, m, n, a, lda, b, ldb);
		}
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void fill(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar val)
{
    launch_fill_kernel(uplo, m, n, a, lda, val);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void scale(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar alpha)
{
	launch_scale_kernel<T_Scalar>(uplo, m, n, a, lda, alpha);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normOne(::cla3p::prop_t ptype, 
                                                 ::cla3p::uplo_t uplo,
                                                 int_t m, int_t n, 
                                                 const T_Scalar* a, 
                                                 int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_one_norm_kernel(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normOne is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normInf(::cla3p::prop_t ptype, 
										         ::cla3p::uplo_t uplo,
										         int_t m, int_t n, 
										         const T_Scalar* a, 
										         int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_inf_norm_kernel(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normInf is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normMax(::cla3p::prop_t ptype, 
	                                             ::cla3p::uplo_t uplo,
                                                 int_t m, int_t n, 
                                                 const T_Scalar* a, 
                                                 int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_max_norm_kernel(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normMax is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normFro(::cla3p::prop_t ptype, 
                                                 ::cla3p::uplo_t uplo,
                                                 int_t m, int_t n, 
                                                 const T_Scalar* a, 
                                                 int_t lda)
{
	typename TypeTraits<T_Scalar>::real_type ret = 0;
	if(ptype == ::cla3p::prop_t::General && uplo == ::cla3p::uplo_t::Full) {
		ret = launch_matrix_fro_norm_kernel(m, n, a, lda);
	} else {
		// TODO: implement
		throw err::CudaException("blk::dns::normFro is only implemented for full matrices.");
	}
	return ret;
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void conjugate(::cla3p::uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda)
{
	launch_conjugate_kernel<T_Scalar>(uplo, m, n, a, lda);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void transpose(int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar* b, int_t ldb)
{
	T_Scalar alpha = makeScalar<T_Scalar>(1);
	T_Scalar beta  = makeScalar<T_Scalar>(0);
    T_Scalar *p_null = nullptr;
	globalCuBlasHandler().geam(::cla3p::op_t::T,
					           ::cla3p::op_t::N,
								n, m,
					            &alpha, a, lda,
					            &beta, p_null, n,
					            b, ldb);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void ctranspose(int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar* b, int_t ldb)
{
	T_Scalar alpha = makeScalar<T_Scalar>(1);
	T_Scalar beta  = makeScalar<T_Scalar>(0);
    T_Scalar *p_null = nullptr;
	globalCuBlasHandler().geam<T_Scalar>(::cla3p::op_t::C,
					                     ::cla3p::op_t::N,
								         n, m,
					                     &alpha, a, lda,
					                     &beta, p_null, n,
					                     b, ldb);
}
/*-------------------------------------------------*/
//
// Update: C += alpha * op(A), C (m x n)
//
template <typename T_Scalar>
void update(int_t m, int_t n, ::cla3p::op_t opA, T_Scalar alpha, const T_Scalar *a, int_t lda, T_Scalar* c, int_t ldc)
{
	T_Scalar beta = makeScalar<T_Scalar>(1);
	globalCuBlasHandler().geam<T_Scalar>(opA,
					                     ::cla3p::op_t::N,
								         m, n,
					                     &alpha, a, lda,
					                     &beta, c, ldc,
					                     c, ldc);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void getReal(::cla3p::uplo_t uplo, 
	           int_t m, int_t n, 
			   const T_Scalar* a, int_t lda, 
			   typename TypeTraits<T_Scalar>::real_type* b, int_t ldb)
{
	launch_get_real_kernel(uplo, m, n, a, lda, b, ldb);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void getImag(::cla3p::uplo_t uplo, 
	            int_t m, int_t n, 
				const T_Scalar* a, int_t lda, 
				typename TypeTraits<T_Scalar>::real_type* b, int_t ldb)
{
	launch_get_imag_kernel(uplo, m, n, a, lda, b, ldb);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void geevCalculateComplexEigenvectors(int_t n, 
                                      const typename TypeTraits<T_Scalar>::real_type* w, 
                                      const typename TypeTraits<T_Scalar>::real_type* vr, int_t ldvr,
                                      T_Scalar *vc, int_t ldvc)
{
    launch_geev_calculate_complex_eigenvectors_kernel(n, w, vr, ldvr, vc, ldvc);
}

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS2D_HPP_
