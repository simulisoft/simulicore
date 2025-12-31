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

#ifndef CLA3P_ALGEBRA_FUNCTIONAL_UPDATE_HPP_
#define CLA3P_ALGEBRA_FUNCTIONAL_UPDATE_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/

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
 * @ingroup cla3p_module_index_math_op_add
 * @brief Update a dense vector with a compatible scaled dense vector.
 * @details Performs the operation <b>Y = Y + alpha * X</b>
 * @param[in] alpha The scaling coefficient.
 * @param[in] X The input dense vector.
 * @param[in,out] Y The dense vector to be updated.
 */
template <typename T_Scalar>
void update(T_Scalar alpha,
    const dns::XxVector<T_Scalar>& X,
    dns::XxVector<T_Scalar>& Y);

/**
 * @ingroup cla3p_module_index_math_op_add
 * @brief Update a dense matrix with a compatible scaled dense matrix.
 * @details Performs the operation <b>B = B + alpha * A</b>
 * @param[in] alpha The scaling coefficient.
 * @param[in] A The input dense matrix.
 * @param[in,out] B The dense matrix to be updated.
 */
template <typename T_Scalar>
void update(T_Scalar alpha,
    const dns::XxMatrix<T_Scalar>& A,
    dns::XxMatrix<T_Scalar>& B);

/**
 * @ingroup cla3p_module_index_math_op_add
 * @brief Update a sparse matrix with a compatible scaled sparse matrix.
 * @details Performs the operation <b>B = B + alpha * A</b>
 * @param[in] alpha The scaling coefficient.
 * @param[in] A The input sparse matrix.
 * @param[in,out] B The sparse matrix to be updated.
 */
template <typename T_Int, typename T_Scalar>
void update(T_Scalar alpha,
    const csc::XxMatrix<T_Int,T_Scalar>& A,
    csc::XxMatrix<T_Int,T_Scalar>& B);

/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_ALGEBRA_FUNCTIONAL_UPDATE_HPP_
