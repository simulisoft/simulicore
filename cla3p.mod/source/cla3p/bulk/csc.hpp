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

#ifndef CLA3P_BULK_CSC_HPP_
#define CLA3P_BULK_CSC_HPP_

/**
 * @file
 */

#include <iostream>

#include "cla3p/types.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace blk {
namespace csc {
/*-------------------------------------------------*/

void roll(int_t n, int_t *colptr);

void unroll(int_t n, int_t *colptr);

int_t maxrlen(int_t n, const int_t *colptr);

void check(prop_t ptype, uplo_t uplo, int_t m, int_t n, const int_t *colptr, const int_t *rowidx);

void sort(int_t n, const int_t *colptr, int_t *rowidx);

template <typename T_Scalar>
void sort(int_t n, const int_t *colptr, int_t *rowidx, T_Scalar *values);

template <typename T_Scalar>
void print_to_stream(std::ostream& os, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values, std::streamsize prec = 0);

template <typename T_Scalar>
void print(int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values, std::streamsize prec = 0)
{
	print_to_stream(std::cout, n, colptr, rowidx, values, prec);
}

template <typename T_Scalar>
void transpose(int_t m, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values,
		int_t *colptr_out, int_t *rowidx_out, T_Scalar *values_out, T_Scalar coeff = 1);

template <typename T_Scalar>
void conjugate_transpose(int_t m, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values,
		int_t *colptr_out, int_t *rowidx_out, T_Scalar *values_out, T_Scalar coeff = 1);

void uplo2ge_colptr(uplo_t uplo, int_t n, const int_t *colptr, const int_t *rowidx, int_t *colptr_out);

template <typename T_Scalar>
void sy2ge(uplo_t uplo, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values,
		int_t *colptr_out, int_t *rowidx_out, T_Scalar *values_out);

template <typename T_Scalar>
void he2ge(uplo_t uplo, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values,
		int_t *colptr_out, int_t *rowidx_out, T_Scalar *values_out);

template <typename T_Scalar>
void remove_duplicates(int_t n, int_t *colptr, int_t *rowidx, T_Scalar *values, dup_t op);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_one(prop_t ptype, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_inf(prop_t ptype, int_t m, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_max(int_t n, const int_t *colptr, const T_Scalar *values);

template <typename T_Scalar>
typename TypeTraits<T_Scalar>::real_type 
norm_fro(prop_t ptype, int_t n, const int_t *colptr, const int_t *rowidx, const T_Scalar *values);

template <typename T_Scalar>
void permute(prop_t ptype, uplo_t uplo, int_t m, int_t n, 
		const int_t *colptr, const int_t *rowidx, const T_Scalar *values,
		int_t *colptr_out, int_t *rowidx_out, T_Scalar *values_out, const int_t *P, const int_t *Q);

/*-------------------------------------------------*/
} // namespace csc
} // namespace blk
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BULK_CSC_HPP_
