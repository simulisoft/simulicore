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

#ifndef CULITE_OPERATORS_ADD_HPP_
#define CULITE_OPERATORS_ADD_HPP_

/**
 * @file
 */

#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"
#include "culite/sparse/csr_xxmatrix.hpp"
#include "culite/sparse/csc_xxmatrix.hpp"

#include "culite/virtuals/virtual_expression.hpp"
#include "culite/virtuals/virtual_object.hpp"
#include "culite/virtuals/virtual_plus_minus.hpp"

/*-------------------------------------------------*/

//
// Virtual + Virtual
//

template <typename T_Result, typename T_Left, typename T_Right>
culite::VirtualPlus<T_Result,T_Left,T_Right>
operator+(const culite::VirtualExpression<T_Result, T_Left>& A, 
	      const culite::VirtualExpression<T_Result, T_Right>& B)
{
	return culite::VirtualPlus<T_Result,T_Left,T_Right>(A, B);
}

/*-------------------------------------------------*/

//
// XxObject + XxObject
//

/**
 * @ingroup culite_module_index_math_operators_add
 * @brief Adds two dense vectors.
 * @details Performs the operation @f$ x + y @f$. Vectors must have compatible sizes.
 * @param[in] x The first vector.
 * @param[in] y The second vector.
 * @return The virtual vector sum representing the addition.
 */
template <typename T_Scalar>
culite::alias::VirtualPlus_vec<T_Scalar>
operator+(const culite::dns::XxVector<T_Scalar>& x,
		  const culite::dns::XxVector<T_Scalar>& y)
{
	return (x.virtualize() + y.virtualize());
}

/**
 * @ingroup culite_module_index_math_operators_add
 * @brief Adds two dense matrices.
 * @details Performs the operation @f$ A + B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first matrix.
 * @param[in] B The second matrix.
 * @return The virtual matrix sum representing the addition.
 */
template <typename T_Scalar>
culite::alias::VirtualPlus_dns<T_Scalar>
operator+(const culite::dns::XxMatrix<T_Scalar>& A,
		  const culite::dns::XxMatrix<T_Scalar>& B)
{
	return (A.virtualize() + B.virtualize());
}

/**
 * @ingroup culite_module_index_math_operators_add
 * @brief Adds two sparse matrices.
 * @details Performs the operation @f$ A + B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first sparse matrix.
 * @param[in] B The second sparse matrix.
 * @return The virtual matrix sum representing the addition.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualPlus_csr<T_Int,T_Scalar>
operator+(const culite::csr::XxMatrix<T_Int,T_Scalar>& A,
		  const culite::csr::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A.virtualize() + B.virtualize());
}

/**
 * @ingroup culite_module_index_math_operators_add
 * @brief Adds two sparse matrices.
 * @details Performs the operation @f$ A + B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first sparse matrix.
 * @param[in] B The second sparse matrix.
 * @return The virtual matrix sum representing the addition.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualPlus_csc<T_Int,T_Scalar>
operator+(const culite::csc::XxMatrix<T_Int,T_Scalar>& A,
		  const culite::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A.virtualize() + B.virtualize());
}

/*-------------------------------------------------*/

//
// XxObject + Virtual
//

template <typename T_Scalar, typename T_Right>
culite::VirtualPlus<
	culite::dns::XxVector<T_Scalar>,
	culite::alias::VirtualObj_vec<T_Scalar>,
	T_Right>
operator+(
		const culite::dns::XxVector<T_Scalar>& A,
		const culite::alias::VirtualExpr_vec<T_Scalar,T_Right>& B)
{
	return (A.virtualize() + B);
}

template <typename T_Scalar, typename T_Right>
culite::VirtualPlus<
	culite::dns::XxMatrix<T_Scalar>,
	culite::alias::VirtualObj_dns<T_Scalar>,
	T_Right>
operator+(
		const culite::dns::XxMatrix<T_Scalar>& A,
		const culite::alias::VirtualExpr_dns<T_Scalar,T_Right>& B)
{
	return (A.virtualize() + B);
}

template <typename T_Int, typename T_Scalar, typename T_Right>
culite::VirtualPlus<
	culite::csr::XxMatrix<T_Int,T_Scalar>,
	culite::alias::VirtualObj_csr<T_Int,T_Scalar>,
	T_Right>
operator+(
		const culite::csr::XxMatrix<T_Int,T_Scalar>& A,
		const culite::alias::VirtualExpr_csr<T_Int,T_Scalar,T_Right>& B)
{
	return (A.virtualize() + B);
}

template <typename T_Int, typename T_Scalar, typename T_Right>
culite::VirtualPlus<
	culite::csc::XxMatrix<T_Int,T_Scalar>,
	culite::alias::VirtualObj_csc<T_Int,T_Scalar>,
	T_Right>
operator+(
		const culite::csc::XxMatrix<T_Int,T_Scalar>& A,
		const culite::alias::VirtualExpr_csc<T_Int,T_Scalar,T_Right>& B)
{
	return (A.virtualize() + B);
}

/*-------------------------------------------------*/

//
// Virtual + XxObject
//

template <typename T_Scalar, typename T_Left>
culite::VirtualPlus<
	culite::dns::XxVector<T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_vec<T_Scalar>>
operator+(
		const culite::alias::VirtualExpr_vec<T_Scalar, T_Left>& A,
		const culite::dns::XxVector<T_Scalar>& B)
{
	return (A + B.virtualize());
}

template <typename T_Scalar, typename T_Left>
culite::VirtualPlus<
	culite::dns::XxMatrix<T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_dns<T_Scalar>>
operator+(
		const culite::alias::VirtualExpr_dns<T_Scalar, T_Left>& A,
		const culite::dns::XxMatrix<T_Scalar>& B)
{
	return (A + B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
culite::VirtualPlus<
	culite::csr::XxMatrix<T_Int,T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_csr<T_Int,T_Scalar>>
operator+(
		const culite::alias::VirtualExpr_csr<T_Int,T_Scalar,T_Left>& A,
		const culite::csr::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A + B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
culite::VirtualPlus<
	culite::csc::XxMatrix<T_Int,T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_csc<T_Int,T_Scalar>>
operator+(
		const culite::alias::VirtualExpr_csc<T_Int,T_Scalar,T_Left>& A,
		const culite::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A + B.virtualize());
}

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

//
// Virtual - Virtual
//

template <typename T_Result, typename T_Left, typename T_Right>
culite::VirtualMinus<T_Result,T_Left,T_Right>
operator-(
	const culite::VirtualExpression<T_Result, T_Left>& A, 
	const culite::VirtualExpression<T_Result, T_Right>& B)
{
	return culite::VirtualMinus<T_Result,T_Left,T_Right>(A, B);
}

/*-------------------------------------------------*/

//
// XxObject - XxObject
//

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Subtracts two dense vectors.
 * @details Performs the operation @f$ x - y @f$. Vectors must have compatible sizes.
 * @param[in] x The first vector.
 * @param[in] y The second vector.
 * @return The virtual vector difference representing the subtraction.
 */
template <typename T_Scalar>
culite::alias::VirtualMinus_vec<T_Scalar>
operator-(
		const culite::dns::XxVector<T_Scalar>& x,
		const culite::dns::XxVector<T_Scalar>& y)
{
	return (x.virtualize() - y.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Subtracts two dense matrices.
 * @details Performs the operation @f$ A - B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first matrix.
 * @param[in] B The second matrix.
 * @return The virtual matrix difference representing the subtraction.
 */
template <typename T_Scalar>
culite::alias::VirtualMinus_dns<T_Scalar>
operator-(
		const culite::dns::XxMatrix<T_Scalar>& A,
		const culite::dns::XxMatrix<T_Scalar>& B)
{
	return (A.virtualize() - B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Subtracts two sparse matrices.
 * @details Performs the operation @f$ A - B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first sparse matrix.
 * @param[in] B The second sparse matrix.
 * @return The virtual matrix difference representing the subtraction.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualMinus_csr<T_Int,T_Scalar>
operator-(
		const culite::csr::XxMatrix<T_Int,T_Scalar>& A,
		const culite::csr::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A.virtualize() - B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Subtracts two sparse matrices.
 * @details Performs the operation @f$ A - B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first sparse matrix.
 * @param[in] B The second sparse matrix.
 * @return The virtual matrix difference representing the subtraction.
 */
template <typename T_Int, typename T_Scalar>
culite::alias::VirtualMinus_csc<T_Int,T_Scalar>
operator-(
		const culite::csc::XxMatrix<T_Int,T_Scalar>& A,
		const culite::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A.virtualize() - B.virtualize());
}

/*-------------------------------------------------*/

//
// XxObject - Virtual
//

template <typename T_Scalar, typename T_Right>
culite::VirtualMinus<
	culite::dns::XxVector<T_Scalar>,
	culite::alias::VirtualObj_vec<T_Scalar>,
	T_Right>
operator-(
		const culite::dns::XxVector<T_Scalar>& A,
		const culite::alias::VirtualExpr_vec<T_Scalar,T_Right>& B)
{
	return (A.virtualize() - B);
}

template <typename T_Scalar, typename T_Right>
culite::VirtualMinus<
	culite::dns::XxMatrix<T_Scalar>,
	culite::alias::VirtualObj_dns<T_Scalar>,
	T_Right>
operator-(
		const culite::dns::XxMatrix<T_Scalar>& A,
		const culite::alias::VirtualExpr_dns<T_Scalar,T_Right>& B)
{
	return (A.virtualize() - B);
}

template <typename T_Int, typename T_Scalar, typename T_Right>
culite::VirtualMinus<
	culite::csr::XxMatrix<T_Int,T_Scalar>,
	culite::alias::VirtualObj_csr<T_Int,T_Scalar>,
	T_Right>
operator-(
		const culite::csr::XxMatrix<T_Int,T_Scalar>& A,
		const culite::alias::VirtualExpr_csr<T_Int,T_Scalar,T_Right>& B)
{
	return (A.virtualize() - B);
}

template <typename T_Int, typename T_Scalar, typename T_Right>
culite::VirtualMinus<
	culite::csc::XxMatrix<T_Int,T_Scalar>,
	culite::alias::VirtualObj_csc<T_Int,T_Scalar>,
	T_Right>
operator-(
		const culite::csc::XxMatrix<T_Int,T_Scalar>& A,
		const culite::alias::VirtualExpr_csc<T_Int,T_Scalar,T_Right>& B)
{
	return (A.virtualize() - B);
}

/*-------------------------------------------------*/

//
// Virtual + XxObject
//

template <typename T_Scalar, typename T_Left>
culite::VirtualMinus<
	culite::dns::XxVector<T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_vec<T_Scalar>>
operator-(
		const culite::alias::VirtualExpr_vec<T_Scalar,T_Left>& A,
		const culite::dns::XxVector<T_Scalar>& B)
{
	return (A - B.virtualize());
}

template <typename T_Scalar, typename T_Left>
culite::VirtualMinus<
	culite::dns::XxMatrix<T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_dns<T_Scalar>>
operator-(
		const culite::alias::VirtualExpr_dns<T_Scalar,T_Left>& A,
		const culite::dns::XxMatrix<T_Scalar>& B)
{
	return (A - B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
culite::VirtualMinus<
	culite::csr::XxMatrix<T_Int,T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_csr<T_Int,T_Scalar>>
operator-(
		const culite::alias::VirtualExpr_csr<T_Int,T_Scalar,T_Left>& A,
		const culite::csr::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A - B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
culite::VirtualMinus<
	culite::csc::XxMatrix<T_Int,T_Scalar>,
	T_Left,
	culite::alias::VirtualObj_csc<T_Int,T_Scalar>>
operator-(
		const culite::alias::VirtualExpr_csc<T_Int,T_Scalar,T_Left>& A,
		const culite::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A - B.virtualize());
}

/*-------------------------------------------------*/

#endif // CULITE_OPERATORS_ADD_HPP_
