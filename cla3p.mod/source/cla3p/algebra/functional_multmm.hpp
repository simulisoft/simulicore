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

#ifndef CLA3P_ALGEBRA_FUNCTIONAL_MULTMM_HPP_
#define CLA3P_ALGEBRA_FUNCTIONAL_MULTMM_HPP_

/**
 * @file
 */

#include "cla3p/types/enums.hpp"
#include "cla3p/sparse/csr_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

namespace cla3p { namespace dns { template <typename T_Scalar> class XxMatrix; } } // namespace cla3p::dns

/*-------------------------------------------------*/
namespace cla3p { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Updates a general matrix with a matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 *          @p C can be General/Symmetric/Hermitian.
 *          If @p C is Symmetric/Hermitian, it is assumed that the outcome of the operation @f$ op_A(A) \cdot op_B(B) @f$ is also Symmetric/Hermitian
 *          and only the corresponding part (upper/lower) will be calculated.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB | C |
 *          |---|---|-----|-----|---|
 *          | General | General | unconstrained | unconstrained | General/Symmetric/Hermitian |
 *          | Symmetric | General | ignored | must be set to N | General |
 *          | Hermitian | General | ignored | must be set to N | General |
 *          | Triangular | General | unconstrained | must be set to N | General |
 *          | General | Symmetric | must be set to N | ignored | General |
 *          | General | Hermitian | must be set to N | ignored | General |
 *          | General | Triangular | must be set to N | unconstrained | General |
 *
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input matrix.
 * @param[in] opB The operation to be performed for matrix @p B.
 * @param[in] B The input matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The matrix to be updated.
 */
template <typename T_Scalar>
void mult(T_Scalar alpha,
          op_t opA, const dns::XxMatrix<T_Scalar>& A,
          op_t opB, const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Replaces a matrix with a scaled triangular matrix-matrix product.
 * @details Performs the operation @f$ B = \alpha \cdot op_A(A) \cdot B @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input triangular matrix.
 * @param[in,out] B The matrix to be replaced.
 */
template <typename T_Scalar>
void trimult(T_Scalar alpha, op_t opA,
             const dns::XxMatrix<T_Scalar>& A,
             dns::XxMatrix<T_Scalar>& B);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Replaces a matrix with a scaled triangular matrix-matrix product.
 * @details Performs the operation @f$ B = \alpha \cdot B \cdot op_A(A) @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in,out] B The matrix to be replaced.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input triangular matrix.
 */
template <typename T_Scalar>
void trimult(T_Scalar alpha,
             dns::XxMatrix<T_Scalar>& B,
             op_t opA, const dns::XxMatrix<T_Scalar>& A);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Replaces a matrix with the scaled solution of a triangular system.
 * @details Solves the system @f$ op_A(A) \cdot X = \alpha \cdot B @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input triangular matrix.
 * @param[in,out] B On entry, the rhs, on exit the system solution @p X.
 */
template <typename T_Scalar>
void trisol(T_Scalar alpha, op_t opA,
            const dns::XxMatrix<T_Scalar>& A,
            dns::XxMatrix<T_Scalar>& B);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Replaces a matrix with the scaled solution of a triangular system.
 * @details Solves the system @f$ X \cdot op_A(A) = \alpha \cdot B @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in,out] B On entry, the rhs, on exit the system solution @p X.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input triangular matrix.
 */
template <typename T_Scalar>
void trisol(T_Scalar alpha,
            dns::XxMatrix<T_Scalar>& B,
            op_t opA, const dns::XxMatrix<T_Scalar>& A);

/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Updates a general dense matrix with a sparse-dense matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot B @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | C |
 *          |---|---|-----|---|
 *          | General | General | unconstrained | General |
 *          | Symmetric | General | ignored | General |
 *          | Hermitian | General | ignored | General |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] B The input dense matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The dense matrix to be updated.
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA, 
          const csr::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Updates a general dense matrix with a sparse-dense matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot B @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | C |
 *          |---|---|-----|---|
 *          | General | General | unconstrained | General |
 *          | Symmetric | General | ignored | General |
 *          | Hermitian | General | ignored | General |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] B The input dense matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The dense matrix to be updated.
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA, 
          const csc::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxMatrix<T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Updates a dense matrix with a sparse-sparse matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB |
 *          |---|---|-----|-----|
 *          | General | General | unconstrained | unconstrained |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] opB The operation to be performed for matrix @p B.
 * @param[in] B The input dense matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The dense matrix to be updated.
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, 
          op_t opA, const csr::XxMatrix<T_Int,T_Scalar>& A,
          op_t opB, const csr::XxMatrix<T_Int,T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Updates a dense matrix with a sparse-sparse matrix-matrix product.
 * @details Performs the operation @f$ C = \beta \cdot C + \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB |
 *          |---|---|-----|-----|
 *          | General | General | unconstrained | unconstrained |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] opB The operation to be performed for matrix @p B.
 * @param[in] B The input dense matrix.
 * @param[in] beta The scaling coefficient for @p C.
 * @param[in,out] C The dense matrix to be updated.
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, 
          op_t opA, const csc::XxMatrix<T_Int,T_Scalar>& A,
          op_t opB, const csc::XxMatrix<T_Int,T_Scalar>& B,
          T_Scalar beta, dns::XxMatrix<T_Scalar>& C);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Creates a general sparse matrix from a sparse-sparse matrix-matrix product.
 * @details Performs the operation @f$ \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB |
 *          |---|---|-----|-----|
 *          | General | General | N or T | N or T |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] opB The operation to be performed for matrix @p B.
 * @param[in] B The input sparse matrix.
 * @return The matrix @f$ \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 */
template <typename T_Int, typename T_Scalar>
csr::XxMatrix<T_Int,T_Scalar> mult(T_Scalar alpha, 
                                   op_t opA, const csr::XxMatrix<T_Int,T_Scalar>& A,
                                   op_t opB, const csr::XxMatrix<T_Int,T_Scalar>& B);

/**
 * @ingroup cla3p_module_index_math_op_matmat
 * @brief Creates a general sparse matrix from a sparse-sparse matrix-matrix product.
 * @details Performs the operation @f$ \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 *          Valid combinations are the following:
 *
 *          | A | B | opA | opB |
 *          |---|---|-----|-----|
 *          | General | General | N or T | N or T |
 *
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A.
 * @param[in] A The input sparse matrix.
 * @param[in] opB The operation to be performed for matrix @p B.
 * @param[in] B The input sparse matrix.
 * @return The matrix @f$ \alpha \cdot op_A(A) \cdot op_B(B) @f$.
 */
template <typename T_Int, typename T_Scalar>
csc::XxMatrix<T_Int,T_Scalar> mult(T_Scalar alpha, 
                                   op_t opA, const csc::XxMatrix<T_Int,T_Scalar>& A,
                                   op_t opB, const csc::XxMatrix<T_Int,T_Scalar>& B);

/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_ALGEBRA_FUNCTIONAL_MULTMM_HPP_
