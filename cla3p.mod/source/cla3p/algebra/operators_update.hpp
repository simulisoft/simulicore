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

#ifndef CLA3P_OPERATORS_UPDATE_HPP_
#define CLA3P_OPERATORS_UPDATE_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/algebra/functional_update.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csr { template <typename T_Int, typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
} // namespace cla3p
/*-------------------------------------------------*/

/*
 * Generic Update Operator
 */
template <typename T_Virtual> 
void operator+=(
		typename T_Virtual::result_type& A,
		const cla3p::VirtualExpression<typename T_Virtual::result_type, T_Virtual>& vB)
{
	vB.accumulateOnExisting(A, 1.0);
}

/*
 * Generic Update Operator
 */
template <typename T_Virtual> 
void operator-=(
		typename T_Virtual::result_type& A,
		const cla3p::VirtualExpression<typename T_Virtual::result_type, T_Virtual>& vB)
{
	vB.accumulateOnExisting(A, -1.0);
}

/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Adds a dense vector to another.
 * @details Performs the in-place addition @f$ y = y + x @f$.
 * @param[in,out] y The vector to be updated.
 * @param[in] x The vector to add.
 */
template <typename T_Scalar>
void operator+=(
		cla3p::dns::XxVector<T_Scalar>& y,
		const cla3p::dns::XxVector<T_Scalar>& x)
{
	cla3p::ops::update(T_Scalar(1), x, y);
}

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Subtracts a dense vector from another.
 * @details Performs the in-place subtraction @f$ y = y - x @f$.
 * @param[in,out] y The vector to be updated.
 * @param[in] x The vector to subtract.
 */
template <typename T_Scalar>
void operator-=(
		cla3p::dns::XxVector<T_Scalar>& y,
		const cla3p::dns::XxVector<T_Scalar>& x)
{
	cla3p::ops::update(T_Scalar(-1), x, y);
}

/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Adds a dense matrix to another.
 * @details Performs the in-place addition @f$ A = A + B @f$. Matrices must
 *          have compatible dimensions and properties.
 * @param[in,out] A The matrix to be updated.
 * @param[in] B The matrix to add.
 */
template <typename T_Scalar>
void operator+=(
		cla3p::dns::XxMatrix<T_Scalar>& A,
		const cla3p::dns::XxMatrix<T_Scalar>& B)
{
	cla3p::ops::update(T_Scalar(1), B, A);
}

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Subtracts a dense matrix from another.
 * @details Performs the in-place subtraction @f$ A = A - B @f$. Matrices must
 *          have compatible dimensions and properties.
 * @param[in,out] A The matrix to be updated.
 * @param[in] B The matrix to subtract.
 */
template <typename T_Scalar>
void operator-=(
		cla3p::dns::XxMatrix<T_Scalar>& A,
		const cla3p::dns::XxMatrix<T_Scalar>& B)
{
	cla3p::ops::update(T_Scalar(-1), B, A);
}

/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Adds a CSR sparse matrix to another.
 * @details Performs the in-place addition @f$ A = A + B @f$. Matrices must
 *          have compatible dimensions and properties.
 * @param[in,out] A The CSR sparse matrix to be updated.
 * @param[in] B The CSR sparse matrix to add.
 */
template <typename T_Int, typename T_Scalar>
void operator+=(
		cla3p::csr::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::csr::XxMatrix<T_Int,T_Scalar>& B)
{
	cla3p::ops::update(T_Scalar(1), B, A);
}

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Subtracts a CSR sparse matrix from another.
 * @details Performs the in-place subtraction @f$ A = A - B @f$. Matrices must
 *          have compatible dimensions and properties.
 * @param[in,out] A The CSR sparse matrix to be updated.
 * @param[in] B The CSR sparse matrix to subtract.
 */
template <typename T_Int, typename T_Scalar>
void operator-=(
		cla3p::csr::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::csr::XxMatrix<T_Int,T_Scalar>& B)
{
	cla3p::ops::update(T_Scalar(-1), B, A);
}

/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Adds a CSC sparse matrix to another.
 * @details Performs the in-place addition @f$ A = A + B @f$. Matrices must
 *          have compatible dimensions and properties.
 * @param[in,out] A The CSC sparse matrix to be updated.
 * @param[in] B The CSC sparse matrix to add.
 */
template <typename T_Int, typename T_Scalar>
void operator+=(
		cla3p::csc::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	cla3p::ops::update(T_Scalar(1), B, A);
}

/**
 * @ingroup cla3p_module_index_math_operators_update
 * @brief Subtracts a CSC sparse matrix from another.
 * @details Performs the in-place subtraction @f$ A = A - B @f$. Matrices must
 *          have compatible dimensions and properties.
 * @param[in,out] A The CSC sparse matrix to be updated.
 * @param[in] B The CSC sparse matrix to subtract.
 */
template <typename T_Int, typename T_Scalar>
void operator-=(
		cla3p::csc::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	cla3p::ops::update(T_Scalar(-1), B, A);
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_UPDATE_HPP_
