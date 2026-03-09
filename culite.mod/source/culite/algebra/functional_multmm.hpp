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

#ifndef CULITE_ALGEBRA_FUNCTIONAL_MULTMM_HPP_
#define CULITE_ALGEBRA_FUNCTIONAL_MULTMM_HPP_

/**
 * @file
 */

#include <cla3p/types/enums.hpp>

#include "culite/generic/cublas_handler.hpp"
#include "culite/generic/cusparse_handler.hpp"
#include "culite/sparse/csr_xxmatrix.hpp"
#include "culite/sparse/csc_xxmatrix.hpp"

namespace culite { namespace dns { template <typename T_Scalar> class XxMatrix; } } // namespace culite::dns

/*-------------------------------------------------*/
namespace culite { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_matmat
 * @brief Updates a general matrix with a matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 *          @p C can be General/Symmetric/Hermitian.
 *          If @p C is Symmetric/Hermitian, it is assumed that the outcome of the operation @f$ op_A(A) \cdot op_B(B) @f$ is also Symmetric/Hermitian
 *          and only the corresponding part (upper/lower) will be calculated.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB | C |
 *          |---|---|-----|-----|---|
 *          | General | General | unconstrained | unconstrained | General |
 *          | Symmetric | General | N | unconstrained | General |
 *          | Hermitian | General | N | unconstrained | General |
 *          | General | General | N or T | T or N | Symmetric |
 *          | General | General | N or C | C or N | Hermitian |
 *
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input matrix.
 * @param[in] opB The operation to be performed for matrix @p B.
 * @param[in] B The input matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The matrix to be updated.
 * @param[in] cuBlasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 */
template <typename T_Scalar>
void mult(T_Scalar alpha,
          ::cla3p::op_t opA, const dns::XxMatrix<T_Scalar>& A,
          ::cla3p::op_t opB, const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C,
          CuBlasHandler& cuBlasHandler = globalCuBlasHandler());

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_matmat
 * @brief Updates a general dense matrix with a sparse-dense matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot B @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB | C |
 *          |---|---|-----|-----|---|
 *          | General | General | unconstrained | N | General |
 *          | Symmetric | General | N | N | General |
 *          | Hermitian | General | N | N | General |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] B The input dense matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The dense matrix to be updated.
 * @param[in] cuSparseHandler The cuSPARSE handler for GPU operations (defaults to global handler).
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA, 
          const csr::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C,
          CuSparseHandler& cuSparseHandler = globalCuSparseHandler());

/**
 * @ingroup culite_module_index_math_op_matmat
 * @brief Updates a general dense matrix with a sparse-dense matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot B @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB | C |
 *          |---|---|-----|-----|---|
 *          | General | General | unconstrained | N | General |
 *          | Symmetric | General | N | N | General |
 *          | Hermitian | General | N | N | General |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] B The input dense matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The dense matrix to be updated.
 * @param[in] cuSparseHandler The cuSPARSE handler for GPU operations (defaults to global handler).
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA, 
          const csc::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C,
          CuSparseHandler& cuSparseHandler = globalCuSparseHandler());

/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_ALGEBRA_FUNCTIONAL_MULTMM_HPP_
