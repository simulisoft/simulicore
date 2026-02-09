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

#ifndef CULITE_ALGEBRA_FUNCTIONAL_ADD_HPP_
#define CULITE_ALGEBRA_FUNCTIONAL_ADD_HPP_

/**
 * @file
 */

#include <cla3p/types/enums.hpp>

#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"
#include "culite/generic/cublas_handler.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_math_op_add
 * @brief Adds two compatible scaled dense vectors.
 * @details Performs the operation @f$ \alpha \cdot x + \beta \cdot y @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient for @p x.
 * @param[in] x The first input dense vector.
 * @param[in] beta The scaling coefficient for @p y.
 * @param[in] y The second input dense vector.
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 * @return The result of the operation @f$ \alpha \cdot x + \beta \cdot y @f$.
 */
template <typename T_Scalar>
dns::XxVector<T_Scalar> add(T_Scalar alpha, const dns::XxVector<T_Scalar>& x,
		                    T_Scalar beta , const dns::XxVector<T_Scalar>& y,
                            CuBlasHandler& cublasHandler = globalCuBlasHandler());

/**
 * @ingroup culite_module_index_math_op_add
 * @brief Adds two compatible scaled dense matrices with optional transpose operations.
 * @details Performs the operation @f$ \alpha \cdot op_A(A) + \beta \cdot op_B(B) @f$,
 *          where @f$ op_A @f$ and @f$ op_B @f$ can be no-transpose, transpose, or conjugate transpose.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] opA The operation to apply to matrix @p A (no-transpose, transpose, or conjugate transpose).
 * @param[in] alpha The scaling coefficient for @p A.
 * @param[in] A The first input dense matrix.
 * @param[in] opB The operation to apply to matrix @p B (no-transpose, transpose, or conjugate transpose).
 * @param[in] beta The scaling coefficient for @p B.
 * @param[in] B The second input dense matrix.
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 * @return The result of the operation @f$ \alpha \cdot op_A(A) + \beta \cdot op_B(B) @f$.
 */
template <typename T_Scalar>
dns::XxMatrix<T_Scalar> add(::cla3p::op_t opA, T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A,
		                    ::cla3p::op_t opB, T_Scalar beta , const dns::XxMatrix<T_Scalar>& B,
                            CuBlasHandler& cublasHandler = globalCuBlasHandler());

/**
 * @ingroup culite_module_index_math_op_add
 * @brief Adds two compatible scaled dense matrices.
 * @details Performs the operation @f$ \alpha \cdot A + \beta \cdot B @f$.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] alpha The scaling coefficient for @p A.
 * @param[in] A The first input dense matrix.
 * @param[in] beta The scaling coefficient for @p B.
 * @param[in] B The second input dense matrix.
 * @param[in] cublasHandler The cuBLAS handler for GPU operations (defaults to global handler).
 * @return The result of the operation @f$ \alpha \cdot A + \beta \cdot B @f$.
 */
template <typename T_Scalar>
dns::XxMatrix<T_Scalar> add(T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A,
		                    T_Scalar beta , const dns::XxMatrix<T_Scalar>& B,
                            CuBlasHandler& cublasHandler = globalCuBlasHandler())
{
    return add(::cla3p::op_t::N, alpha, A,
               ::cla3p::op_t::N, beta, B,
               cublasHandler);
}


/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_ALGEBRA_FUNCTIONAL_ADD_HPP_
