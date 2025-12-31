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

#ifndef CLA3P_BULK_DNS_HPP_
#define CLA3P_BULK_DNS_HPP_

/**
 * @file
 */

#include "cla3p/types.hpp"

/*-------------------------------------------------*/
namespace cla3p {
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

//
// Get value (column-major)
//
template <typename T_Scalar>
inline T_Scalar& entry(int_t lda, T_Scalar *a, int_t i, int_t j)
{
	return *ptrmv(lda,a,i,j);
}

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

/*-------------------------------------------------*/
} // namespace dns
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BULK_DNS_HPP_
