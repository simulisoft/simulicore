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

#ifndef CULITE_ALGEBRA_FUNCTIONAL_OUTER_HPP_
#define CULITE_ALGEBRA_FUNCTIONAL_OUTER_HPP_

/**
 * @file
 */

#include "culite/generic/cublas_handler.hpp"

namespace culite {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
} // namespace culite

/*-------------------------------------------------*/
namespace culite { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_vecvec
 * @brief Update a matrix with a vector outer product.
 * @details Performs the operation @f$ A = A + \alpha \cdot x y^T @f$.
 *          @p A can be General/Symmetric/Hermitian.
 *          If @p A is Symmetric/Hermitian, it is assumed that the outcome of the operation @f$ x y^T @f$ is also Symmetric/Hermitian
 *          and only the corresponding part (upper/lower) will be calculated.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The update coefficient.
 * @param[in] x The input vector.
 * @param[in] y The input vector.
 * @param[in,out] A The matrix to be updated.
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 */
template <typename T_Scalar>
void outer(T_Scalar alpha, 
           const dns::XxVector<T_Scalar>& x, 
           const dns::XxVector<T_Scalar>& y, 
           dns::XxMatrix<T_Scalar>& A,
           CuBlasHandler& cublasHandler = globalCuBlasHandler());

/**
 * @ingroup culite_module_index_math_op_vecvec
 * @brief Update a matrix with a conjugate vector outer product.
 * @details Performs the operation @f$ A = A + \alpha \cdot x y^H @f$.
 *          @p A can be General/Symmetric/Hermitian.
 *          If @p A is Symmetric/Hermitian, it is assumed that the outcome of the operation @f$ x y^H @f$ is also Symmetric/Hermitian
 *          and only the corresponding part (upper/lower) will be calculated.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The update coefficient.
 * @param[in] x The input vector.
 * @param[in] y The input vector.
 * @param[in,out] A The matrix to be updated.
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 */
template <typename T_Scalar>
void outerc(T_Scalar alpha,
            const dns::XxVector<T_Scalar>& x,
            const dns::XxVector<T_Scalar>& y,
            dns::XxMatrix<T_Scalar>& A,
            CuBlasHandler& cublasHandler = globalCuBlasHandler());

/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_ALGEBRA_FUNCTIONAL_OUTER_HPP_
