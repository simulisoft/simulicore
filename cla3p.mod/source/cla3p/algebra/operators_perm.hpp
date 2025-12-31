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

#ifndef CLA3P_OPERATORS_PERM_HPP_
#define CLA3P_OPERATORS_PERM_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/perms/pxmatrix.hpp"

/*-------------------------------------------------*/
//namespace cla3p {
//namespace prm { template <typename T_Int> class PxMatrix; }
//namespace dns { template <typename T_Scalar> class XxVector; }
//namespace dns { template <typename T_Scalar> class XxMatrix; }
//namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
//} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Multiplies permutation matrix with a vector.
 * @details Performs the operation <b>P * X</b>
 * @param[in] P The input permutation matrix.
 * @param[in] X The input vector.
 * @return The permuted vector.
 */
template <typename T_Int, typename T_Scalar>
cla3p::dns::XxVector<T_Scalar> operator*(
	const cla3p::prm::PxMatrix<T_Int>& P, 
	const cla3p::dns::XxVector<T_Scalar>& X);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Multiplies permutation matrix with a general dense matrix.
 * @details Performs the operation <b>P * A</b>
 * @param[in] P The input permutation matrix.
 * @param[in] A The input general matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::dns::XxMatrix<T_Scalar> operator*(
	const cla3p::prm::PxMatrix<T_Int>& P, 
	const cla3p::dns::XxMatrix<T_Scalar>& A);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Multiplies general dense matrix with a permutation matrix.
 * @details Performs the operation <b>A * P</b>
 * @param[in] A The input general matrix.
 * @param[in] P The input permutation matrix.
 * @return The permuted matrix.
 */

template <typename T_Int, typename T_Scalar>
cla3p::dns::XxMatrix<T_Scalar> operator*(
	const cla3p::dns::XxMatrix<T_Scalar>& A, 
	const cla3p::prm::PxMatrix<T_Int>& P);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Multiplies permutation matrix with a general sparse matrix.
 * @details Performs the operation <b>P * A</b>
 * @param[in] P The input permutation matrix.
 * @param[in] A The input general matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::csc::XxMatrix<T_Int,T_Scalar> operator*(
	const cla3p::prm::PxMatrix<T_Int>& P, 
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Multiplies general sparse matrix with a permutation matrix.
 * @details Performs the operation <b>A * P</b>
 * @param[in] A The input general matrix.
 * @param[in] P The input permutation matrix.
 * @return The permuted matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::csc::XxMatrix<T_Int,T_Scalar> operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::prm::PxMatrix<T_Int>& P);

/**
 * @ingroup cla3p_module_index_math_operators_perm
 * @brief Multiplies a permutation matrix with a permutation matrix.
 * @details Performs the operation <b>P * Q</b>
 * @param[in] P The input permutation matrix.
 * @param[in] Q The input permutation matrix.
 * @return The permuted permutation matrix.
 */
template <typename T_Int>
cla3p::prm::PxMatrix<T_Int> operator*(
	const cla3p::prm::PxMatrix<T_Int>& P,
	const cla3p::prm::PxMatrix<T_Int>& Q);

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

template <typename T_Int, typename T_Virtual>
typename T_Virtual::result_type operator*(
	const cla3p::prm::PxMatrix<T_Int>& P, 
	const cla3p::VirtualExpression<typename T_Virtual::result_type,T_Virtual>& v)
{
	return (P * v.evaluate());
}

template <typename T_Int, typename T_Virtual>
typename T_Virtual::result_type operator*(
	const cla3p::VirtualExpression<typename T_Virtual::result_type,T_Virtual>& v,
	const cla3p::prm::PxMatrix<T_Int>& P) 
{
	return (v.evaluate() * P);
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_PERM_HPP_
