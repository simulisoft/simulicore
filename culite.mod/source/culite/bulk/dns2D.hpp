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
	if(uplo == ::cla3p::uplo_t::Full && lda == m) {
		globalCuBlasHandler().scal<T_Scalar>(m*n, &coeff, a, 1);
	} else {
		// TODO: implement
		err::CudaException("blk::dns::scale2D is not implemented yet for discontiguous memory.");
	}
}
/*-------------------------------------------------*/
//
// Norm 1
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normOne2D(
	::cla3p::prop_t /*ptype*/, 
	::cla3p::uplo_t /*uplo*/,
	int_t /*m*/, int_t /*n*/, 
	const T_Scalar* /*a*/, 
	int_t /*lda*/)
{
	// TODO: implement
	err::CudaException("blk::dns::normOne2D is not implemented yet.");
	return typename TypeTraits<T_Scalar>::real_type(0);
}
/*-------------------------------------------------*/
//
// Norm Inf
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normInf2D(
	::cla3p::prop_t /*ptype*/, 
	::cla3p::uplo_t /*uplo*/,
	int_t /*m*/, int_t /*n*/, 
	const T_Scalar* /*a*/, 
	int_t /*lda*/)
{
	// TODO: implement
	err::CudaException("blk::dns::normInf2D is not implemented yet.");
	return typename TypeTraits<T_Scalar>::real_type(0);
}
/*-------------------------------------------------*/
//
// Norm Max
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normMax2D(
	::cla3p::prop_t /*ptype*/, 
	::cla3p::uplo_t /*uplo*/,
	int_t /*m*/, int_t /*n*/, 
	const T_Scalar* /*a*/, 
	int_t /*lda*/)
{
	// TODO: implement
	err::CudaException("blk::dns::normMax2D is not implemented yet.");
	return typename TypeTraits<T_Scalar>::real_type(0);
}
/*-------------------------------------------------*/
//
// Norm Frobenius
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type normFro2D(
	::cla3p::prop_t /*ptype*/, 
	::cla3p::uplo_t /*uplo*/,
	int_t /*m*/, int_t /*n*/, 
	const T_Scalar* /*a*/, 
	int_t /*lda*/)
{
	// TODO: implement
	err::CudaException("blk::dns::normFro2D is not implemented yet.");
	return typename TypeTraits<T_Scalar>::real_type(0);
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
	err::CudaException("blk::dns::conjugate2D is not implemented yet.");
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
	err::CudaException("blk::dns::getReal2D is not implemented yet.");
}
/*-------------------------------------------------*/
//
// get Imag 2D
//
template <typename T_Scalar>
void getImag2D(::cla3p::uplo_t, int_t, int_t, const T_Scalar*, int_t, typename TypeTraits<T_Scalar>::real_type*, int_t)
{
	// TODO: implement
	err::CudaException("blk::dns::getImag2D is not implemented yet.");
}
/*-------------------------------------------------*/


















/*-------------------------------------------------*/

#if 0

//
// Set zeros on diagonal depending on property
//
template <typename T_Scalar>
void set_diag_zeros(prop_t ptype, int_t n, T_Scalar *a, int_t lda);

//
// Set all entries to val
//
template <typename T_Scalar>
void fill(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar val);

//
// Set all entries to val (separate diagonal)
//
template <typename T_Scalar>
void fill(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar val, T_Scalar dval);

//
// Set all entries to zero
//
template <typename T_Scalar>
void zero(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda)
{
	fill(uplo, m, n, a, lda, T_Scalar(0));
}

//
// Set random values
//
template <typename T_Scalar>
void rand(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda,
		typename TypeTraits<T_Scalar>::real_type lo,
		typename TypeTraits<T_Scalar>::real_type hi);

//
// Copy
//
template <typename T_Scalar>
void copy(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda,
		T_Scalar *b, int_t ldb, T_Scalar coeff = T_Scalar(1));

//
// Get real part from complex
//
template <typename T_Scalar>
void get_real(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda,
		typename TypeTraits<T_Scalar>::real_type *b, int_t ldb);

//
// Get imag part from complex
//
template <typename T_Scalar>
void get_imag(uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda,
		typename TypeTraits<T_Scalar>::real_type *b, int_t ldb);

//
// Set real part to complex
//
template <typename T_Scalar>
void set_real(uplo_t uplo, int_t m, int_t n, const typename TypeTraits<T_Scalar>::real_type *a, int_t lda, 
		T_Scalar *b, int_t ldb);

//
// Set imag part to complex
//
template <typename T_Scalar>
void set_imag(uplo_t uplo, int_t m, int_t n, const typename TypeTraits<T_Scalar>::real_type *a, int_t lda,
		T_Scalar *b, int_t ldb);

//
// Scale
//
template <typename T_Scalar>
void scale(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar coeff);

//
// Transpositions
//
template <typename T_Scalar>
void transpose(int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar *b, int_t ldb, T_Scalar coeff = T_Scalar(1));

//
// Conjugate transpositions
//
template <typename T_Scalar>
void conjugate_transpose(int_t m, int_t n, const T_Scalar *a, int_t lda, T_Scalar *b, int_t ldb, T_Scalar coeff = T_Scalar(1));

//
// Conjugations
//
template <typename T_Scalar>
void conjugate(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda, T_Scalar coeff = T_Scalar(1));

//
// Symmetric to general
//
template <typename T_Scalar>
void sy2ge(uplo_t uplo, int_t n, T_Scalar *a, int_t lda);

//
// Hermitian to general
//
template <typename T_Scalar>
void he2ge(uplo_t uplo, int_t n, T_Scalar *a, int_t lda);

//
// Skew to general
//
template <typename T_Scalar>
void sk2ge(uplo_t uplo, int_t n, T_Scalar *a, int_t lda);

//
// Triangular to general
//
template <typename T_Scalar>
void tr2ge(uplo_t uplo, int_t m, int_t n, T_Scalar *a, int_t lda);

//
// Norm 1
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type norm_one(prop_t ptype, uplo_t uplo,
		int_t m, int_t n, const T_Scalar *a, int_t lda);

//
// Norm Inf
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type norm_inf(prop_t ptype, uplo_t uplo,
		int_t m, int_t n, const T_Scalar *a, int_t lda);

//
// Norm Max
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type norm_max(prop_t ptype, uplo_t uplo,
		int_t m, int_t n, const T_Scalar *a, int_t lda);

//
// Norm Frobenius
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type norm_fro(prop_t ptype, uplo_t uplo,
		int_t m, int_t n, const T_Scalar *a, int_t lda);

//
// Norm Euclidean
//
template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type norm_euc(int_t n, const T_Scalar *a);

//
// Permutations
//
// prop: General                  B = P*A*Q     if P,Q is nullptr, the identity perm is used
// prop: Symmetric/Hermitian/Skew B = P*A*P^{T} Q is not referenced
//
template <typename T_Scalar>
void permute(prop_t ptype, uplo_t uplo, int_t m, int_t n, const T_Scalar *a, int_t lda, 
		T_Scalar *b, int_t ldb, const int_t *P, const int_t *Q);

#endif // 0

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_BULK_DNS_HPP_
