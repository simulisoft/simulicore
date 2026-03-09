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

#ifndef CULITE_ALGEBRA_FUNCTIONAL_MULTMV_HPP_
#define CULITE_ALGEBRA_FUNCTIONAL_MULTMV_HPP_

/**
 * @file
 */

#include <cla3p/types/enums.hpp>

#include "culite/generic/cublas_handler.hpp"
#include "culite/generic/cusparse_handler.hpp"

namespace culite {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csr { template <typename T_Int, typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
} // namespace culite

/*-------------------------------------------------*/
namespace culite { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_matvec
 * @brief Updates a vector with a matrix-vector product.
 * @details Performs the operation @f$ y = \beta \cdot y + \alpha \cdot op_A(A) \cdot x @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A. If @p A is symmetric or hermitian, @p opA is ignored.
 * @param[in] A The input matrix.
 * @param[in] x The input vector.
 * @param[in] beta The scaling coefficient for @p y.
 * @param[in,out] y The vector to be updated.
 * @param[in] cuBlasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 */
template <typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA,
          const dns::XxMatrix<T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta, dns::XxVector<T_Scalar>& y,
          CuBlasHandler& cuBlasHandler = globalCuBlasHandler());

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_matvec
 * @brief Updates a vector with a matrix-vector product.
 * @details Performs the operation @f$ y = \beta \cdot y + \alpha \cdot op_A(A) \cdot x @f$.
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A. If @p A is symmetric or hermitian, @p opA is ignored.
 * @param[in] A The input matrix.
 * @param[in] x The input vector.
 * @param[in] beta The scaling coefficient for @p y.
 * @param[in,out] y The vector to be updated.
 * @param[in] cuSparseHandler The cuSPARSE handler for GPU operations (defaults to global handler).
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA,
          const csr::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta,
          dns::XxVector<T_Scalar>& y,
          CuSparseHandler& cuSparseHandler = globalCuSparseHandler());

/**
 * @ingroup culite_module_index_math_op_matvec
 * @brief Updates a vector with a matrix-vector product.
 * @details Performs the operation @f$ y = \beta \cdot y + \alpha \cdot op_A(A) \cdot x @f$.
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix @p A. If @p A is symmetric or hermitian, @p opA is ignored.
 * @param[in] A The input matrix.
 * @param[in] x The input vector.
 * @param[in] beta The scaling coefficient for @p y.
 * @param[in,out] y The vector to be updated.
 * @param[in] cuSparseHandler The cuSPARSE handler for GPU operations (defaults to global handler).
 */
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, ::cla3p::op_t opA,
          const csc::XxMatrix<T_Int,T_Scalar>& A,
          const dns::XxVector<T_Scalar>& x,
          T_Scalar beta,
          dns::XxVector<T_Scalar>& y,
          CuSparseHandler& cuSparseHandler = globalCuSparseHandler());

/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_ALGEBRA_FUNCTIONAL_MULTMV_HPP_
