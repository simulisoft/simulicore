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

#ifndef CLA3P_ALGEBRA_FUNCTIONAL_ADD_HPP_
#define CLA3P_ALGEBRA_FUNCTIONAL_ADD_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/lra/lra_xxmatrix.hpp"

namespace cla3p{ namespace lra { template <typename T_Matrix> class RankModerator; } }

/*-------------------------------------------------*/
namespace cla3p { 
namespace ops {
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_op_add
 * @brief Adds two compatible scaled dense vectors.
 * @details Performs the operation <b>alpha * X + beta * Y</b>
 * @param[in] alpha The scaling coefficient for X.
 * @param[in] X The first input dense vector.
 * @param[in] beta The scaling coefficient for Y.
 * @param[in] Y The second input dense vector.
 * @return The result of the operation <b>(alpha * X + beta * Y)</b>.
 */
template <typename T_Scalar>
dns::XxVector<T_Scalar> add(
		T_Scalar alpha, const dns::XxVector<T_Scalar>& X,
		T_Scalar beta , const dns::XxVector<T_Scalar>& Y);

/**
 * @ingroup cla3p_module_index_math_op_add
 * @brief Adds two compatible scaled dense matrices.
 * @details Performs the operation <b>alpha * A + beta * B</b>
 * @param[in] alpha The scaling coefficient for A.
 * @param[in] A The first input dense matrix.
 * @param[in] beta The scaling coefficient for B.
 * @param[in] B The second input dense matrix.
 * @return The result of the operation <b>(alpha * A + beta * B)</b>.
 */
template <typename T_Scalar>
dns::XxMatrix<T_Scalar> add(
		T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A,
		T_Scalar beta , const dns::XxMatrix<T_Scalar>& B);

/**
 * @ingroup cla3p_module_index_math_op_add
 * @brief Adds two compatible sparse matrices.
 * @details Performs the operation <b>alpha * A + beta * B</b>
 * @param[in] alpha The scaling coefficient for A.
 * @param[in] A The first input sparse matrix.
 * @param[in] beta The scaling coefficient for B.
 * @param[in] B The second input sparse matrix.
 * @return The result of the operation <b>(alpha * A + beta * B)</b>.
 */
template <typename T_Int, typename T_Scalar>
csc::XxMatrix<T_Int, T_Scalar> add(
		T_Scalar alpha, const csc::XxMatrix<T_Int, T_Scalar>& A,
		T_Scalar beta , const csc::XxMatrix<T_Int, T_Scalar>& B);

/**
 * @ingroup cla3p_module_index_math_op_add
 * @brief Adds two compatible low-rank matrices.
 * @details Performs the operation <b>alpha * A + beta * B</b>
 * @param[in] alpha The scaling coefficient for A.
 * @param[in] A The first input low-rank matrix.
 * @param[in] beta The scaling coefficient for B.
 * @param[in] B The second input low-rank matrix.
 * @param[in] rmod The rank moderator that controls the rank cutoff.
 * @return The result of the operation <b>(alpha * A + beta * B)</b>.
 */
template <typename T_Matrix>
lra::XxMatrix<T_Matrix> add(
		typename T_Matrix::value_type alpha, const lra::XxMatrix<T_Matrix>& A,
		typename T_Matrix::value_type beta, const lra::XxMatrix<T_Matrix>& B,
		const lra::RankModerator<T_Matrix>& rmod);

/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_ALGEBRA_FUNCTIONAL_ADD_HPP_
