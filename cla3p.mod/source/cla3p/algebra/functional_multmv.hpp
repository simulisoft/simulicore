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

#ifndef CLA3P_ALGEBRA_FUNCTIONAL_MULTMV_HPP_
#define CLA3P_ALGEBRA_FUNCTIONAL_MULTMV_HPP_

/**
 * @file
 */

#include "cla3p/types/enums.hpp"

namespace cla3p {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
} // namespace cla3p

/*-------------------------------------------------*/
namespace cla3p { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_op_matvec
 * @brief Updates a vector with a matrix-vector product.
 * @details Performs the operation <b>Y := beta * Y + alpha * opA(A) * X</b>
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix A. If A is symmetric or hermitian, opA is ignored.
 * @param[in] A The input matrix.
 * @param[in] X The input vector.
 * @param[in] beta The scaling coefficient for Y.
 * @param[in,out] Y The vector to be updated.
 */

template <typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
    const dns::XxMatrix<T_Scalar>& A,
    const dns::XxVector<T_Scalar>& X,
		T_Scalar beta,
    dns::XxVector<T_Scalar>& Y);

/**
 * @ingroup cla3p_module_index_math_op_matvec
 * @brief Replaces a vector with a triangular matrix-vector product.
 * @details Performs the operation <b>X := opA(A) * X</b>
 * @param[in] opA The operation to be performed for matrix A.
 * @param[in] A The input triangular matrix.
 * @param[in,out] X The vector to be replaced.
 */
template <typename T_Scalar>
void trimult(op_t opA,
    const dns::XxMatrix<T_Scalar>& A,
    dns::XxVector<T_Scalar>& X);

/**
 * @ingroup cla3p_module_index_math_op_matvec
 * @brief Replaces a vector with the solution of a triangular system.
 * @details Solves the system <b>opA(A) * X = B</b>
 * @param[in] opA The operation to be performed for matrix A.
 * @param[in] A The input triangular matrix.
 * @param[in,out] B On entry, the rhs, on exit the system solution X.
 */
template <typename T_Scalar>
void trisol(op_t opA,
    const dns::XxMatrix<T_Scalar>& A,
    dns::XxVector<T_Scalar>& B);

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_op_matvec
 * @brief Updates a vector with a matrix-vector product.
 * @details Performs the operation <b>Y := beta * Y + alpha * opA(A) * X</b>
 * @param[in] alpha The scaling coefficient.
 * @param[in] opA The operation to be performed for matrix A. If A is symmetric or hermitian, opA is ignored.
 * @param[in] A The input matrix.
 * @param[in] X The input vector.
 * @param[in] beta The scaling coefficient for Y.
 * @param[in,out] Y The vector to be updated.
 */

template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
    const csc::XxMatrix<T_Int,T_Scalar>& A,
    const dns::XxVector<T_Scalar>& X,
		T_Scalar beta,
    dns::XxVector<T_Scalar>& Y);

/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_ALGEBRA_FUNCTIONAL_MULTMV_HPP_
