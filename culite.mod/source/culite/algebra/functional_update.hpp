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

#ifndef CULITE_ALGEBRA_FUNCTIONAL_UPDATE_HPP_
#define CULITE_ALGEBRA_FUNCTIONAL_UPDATE_HPP_

/**
 * @file
 */

#include <cla3p/types/enums.hpp>

#include "culite/generic/cublas_handler.hpp"

/*-------------------------------------------------*/

namespace culite {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
} // namespace culite

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_add
 * @brief Update a dense vector with a compatible scaled dense vector.
 * @details Performs the operation @f$ y = y + \alpha \cdot x @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] x The input dense vector.
 * @param[in,out] y The dense vector to be updated.
 */
template <typename T_Scalar>
void update(T_Scalar alpha,
            const dns::XxVector<T_Scalar>& x,
            dns::XxVector<T_Scalar>& y,
            CuBlasHandler& cublasHandler = globalCuBlasHandler());

template <typename T_Scalar>
void update(::cla3p::op_t opA,
            T_Scalar alpha,
            const dns::XxMatrix<T_Scalar>& A,
            dns::XxMatrix<T_Scalar>& B,
            CuBlasHandler& cublasHandler = globalCuBlasHandler());

/**
 * @ingroup culite_module_index_math_op_add
 * @brief Update a dense matrix with a compatible scaled dense matrix.
 * @details Performs the operation @f$ B = B + \alpha \cdot A @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient.
 * @param[in] A The input dense matrix.
 * @param[in,out] B The dense matrix to be updated.
 */
template <typename T_Scalar>
void update(T_Scalar alpha,
            const dns::XxMatrix<T_Scalar>& A,
            dns::XxMatrix<T_Scalar>& B,
            CuBlasHandler& cublasHandler = globalCuBlasHandler())
{
    update(::cla3p::op_t::N, alpha, A, B);
}

/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_ALGEBRA_FUNCTIONAL_UPDATE_HPP_
