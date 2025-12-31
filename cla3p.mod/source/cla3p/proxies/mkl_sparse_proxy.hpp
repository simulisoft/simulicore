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

#ifndef CLA3P_MKL_SPARSE_PROXY_HPP_
#define CLA3P_MKL_SPARSE_PROXY_HPP_

/**
 * @file
 */

#include "cla3p/types.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace mkl {
/*-------------------------------------------------*/

// C(m x n)
template <typename T_Scalar>
void csc_add(int_t m, int_t n, T_Scalar alpha, op_t opA,
		const int_t* colptrA, const int_t* rowidxA, const T_Scalar* valuesA, 
		const int_t* colptrB, const int_t* rowidxB, const T_Scalar* valuesB,
		int_t **colptrC, int_t **rowidxC, T_Scalar **valuesC);

// A(m x n)
template <typename T_Scalar>
void csc_mv(prop_t propA, uplo_t uploA, int_t m, int_t n, T_Scalar alpha, op_t opA,
		const int_t* colptrA, const int_t* rowidxA, const T_Scalar* valuesA, 
		const T_Scalar* x, T_Scalar beta, T_Scalar *y);

// A(m x n) B(? x k) C(? x k)
template <typename T_Scalar>
void csc_mm(prop_t propA, uplo_t uploA, int_t m, int_t n, T_Scalar alpha, op_t opA,
		const int_t* colptrA, const int_t* rowidxA, const T_Scalar* valuesA, 
		int_t k, const T_Scalar* b, int_t ldb, T_Scalar beta, T_Scalar *c, int_t ldc);

#if 0
template <typename T_Scalar>
void csc_spmm(op_t opA,
    int_t mA, int_t nA, const int_t* colptrA, const int_t* rowidxA, const T_Scalar* valuesA,
    int_t mB, int_t nB, const int_t* colptrB, const int_t* rowidxB, const T_Scalar* valuesB,
    int_t** colptrC, int_t** rowidxC, T_Scalar** valuesC);
#endif

// A(mA x nA) B(mB x nB) C(? x ?)
template <typename T_Scalar>
void csc_sp2m(
    op_t opA, int_t mA, int_t nA, const int_t* colptrA, const int_t* rowidxA, const T_Scalar* valuesA,
    op_t opB, int_t mB, int_t nB, const int_t* colptrB, const int_t* rowidxB, const T_Scalar* valuesB,
    int_t** colptrC, int_t** rowidxC, T_Scalar** valuesC);

// A(mA x nA) B(mB x nB) C(? x ?)
template <typename T_Scalar>
void csc_sp2md(T_Scalar alpha,
    op_t opA, int_t mA, int_t nA, const int_t* colptrA, const int_t* rowidxA, const T_Scalar* valuesA,
    op_t opB, int_t mB, int_t nB, const int_t* colptrB, const int_t* rowidxB, const T_Scalar* valuesB,
    T_Scalar beta, T_Scalar* c, int_t ldc);

/*-------------------------------------------------*/
} // namespace mkl
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_MKL_SPARSE_PROXY_HPP_
