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

#ifndef CULITE_ALGEBRA_FUNCTIONAL_INNER_HPP_
#define CULITE_ALGEBRA_FUNCTIONAL_INNER_HPP_

/**
 * @file
 */

#include <culite/generic/cublas_handler.hpp>

namespace culite { namespace dns { template <typename T_Scalar> class XxVector; } }

/*-------------------------------------------------*/
namespace culite { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_vecvec
 * @brief Calculate vector dot product.
 * @details Performs the operation @f$ x^T \cdot y @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] x The input vector (must have device data).
 * @param[in] y The input vector (must have device data).
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 * @return The dot product value on host.
 */
template <typename T_Scalar>
T_Scalar dot(const dns::XxVector<T_Scalar>& x, const dns::XxVector<T_Scalar>& y, CuBlasHandler& cublasHandler = globalCuBlasHandler());

/**
 * @ingroup culite_module_index_math_op_vecvec
 * @brief Calculate vector conjugate dot product.
 * @details Performs the operation @f$ x^H \cdot y @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] x The input vector (must have device data).
 * @param[in] y The input vector (must have device data).
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 * @return The conjugate dot product value on host.
 */
template <typename T_Scalar>
T_Scalar dotc(const dns::XxVector<T_Scalar>& x, const dns::XxVector<T_Scalar>& y, CuBlasHandler& cublasHandler = globalCuBlasHandler());

/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_ALGEBRA_FUNCTIONAL_INNER_HPP_
