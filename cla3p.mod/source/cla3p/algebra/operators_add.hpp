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

#ifndef CLA3P_OPERATORS_ADD_HPP_
#define CLA3P_OPERATORS_ADD_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_plus_minus.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace alias { 

template <typename T_Scalar>
using VirtualPlus_vec = VirtualPlus<
	dns::XxVector<T_Scalar>,
	VirtualObject<dns::XxVector<T_Scalar>>,
	VirtualObject<dns::XxVector<T_Scalar>>>;

template <typename T_Scalar>
using VirtualPlus_dns = VirtualPlus<
	dns::XxMatrix<T_Scalar>,
	VirtualObject<dns::XxMatrix<T_Scalar>>,
	VirtualObject<dns::XxMatrix<T_Scalar>>>;

template <typename T_Int, typename T_Scalar>
using VirtualPlus_csc = VirtualPlus<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>>;

template <typename T_Scalar>
using VirtualMinus_vec = VirtualMinus<
	dns::XxVector<T_Scalar>,
	VirtualObject<dns::XxVector<T_Scalar>>,
	VirtualObject<dns::XxVector<T_Scalar>>>;

template <typename T_Scalar>
using VirtualMinus_dns = VirtualMinus<
	dns::XxMatrix<T_Scalar>,
	VirtualObject<dns::XxMatrix<T_Scalar>>,
	VirtualObject<dns::XxMatrix<T_Scalar>>>;

template <typename T_Int, typename T_Scalar>
using VirtualMinus_csc = VirtualMinus<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>>;

} // namespace alias
} // namespace cla3p
/*-------------------------------------------------*/

//
// Virtual + Virtual
//

template <typename T_Result, typename T_Left, typename T_Right>
cla3p::VirtualPlus<T_Result,T_Left,T_Right>
operator+(
	const cla3p::VirtualExpression<T_Result, T_Left>& A, 
	const cla3p::VirtualExpression<T_Result, T_Right>& B)
{
	return cla3p::VirtualPlus<T_Result,T_Left,T_Right>(A, B);
}

/*-------------------------------------------------*/

//
// XxObject + XxObject
//

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Adds two dense vectors.
 * @details Performs the operation @f$ x + y @f$. Vectors must have compatible sizes.
 * @param[in] x The first vector.
 * @param[in] y The second vector.
 * @return The virtual vector sum representing the addition.
 */
template <typename T_Scalar>
cla3p::alias::VirtualPlus_vec<T_Scalar>
operator+(
		const cla3p::dns::XxVector<T_Scalar>& x,
		const cla3p::dns::XxVector<T_Scalar>& y)
{
	return (x.virtualize() + y.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Adds two dense matrices.
 * @details Performs the operation @f$ A + B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first matrix.
 * @param[in] B The second matrix.
 * @return The virtual matrix sum representing the addition.
 */
template <typename T_Scalar>
cla3p::alias::VirtualPlus_dns<T_Scalar>
operator+(
		const cla3p::dns::XxMatrix<T_Scalar>& A,
		const cla3p::dns::XxMatrix<T_Scalar>& B)
{
	return (A.virtualize() + B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_add
 * @brief Adds two sparse matrices.
 * @details Performs the operation @f$ A + B @f$. Matrices must have compatible
 *          dimensions and properties.
 * @param[in] A The first sparse matrix.
 * @param[in] B The second sparse matrix.
 * @return The virtual matrix sum representing the addition.
 */
template <typename T_Int, typename T_Scalar>
cla3p::alias::VirtualPlus_csc<T_Int,T_Scalar>
operator+(
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A.virtualize() + B.virtualize());
}

/*-------------------------------------------------*/

//
// XxObject + Virtual
//

template <typename T_Scalar, typename T_Right>
cla3p::VirtualPlus<
	cla3p::dns::XxVector<T_Scalar>,
	cla3p::VirtualObject<cla3p::dns::XxVector<T_Scalar>>,
	T_Right>
operator+(
		const cla3p::dns::XxVector<T_Scalar>& A,
		const cla3p::VirtualExpression<cla3p::dns::XxVector<T_Scalar>, T_Right>& B)
{
	return (A.virtualize() + B);
}

template <typename T_Scalar, typename T_Right>
cla3p::VirtualPlus<
	cla3p::dns::XxMatrix<T_Scalar>,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>,
	T_Right>
operator+(
		const cla3p::dns::XxMatrix<T_Scalar>& A,
		const cla3p::VirtualExpression<cla3p::dns::XxMatrix<T_Scalar>, T_Right>& B)
{
	return (A.virtualize() + B);
}

template <typename T_Int, typename T_Scalar, typename T_Right>
cla3p::VirtualPlus<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>,
	T_Right>
operator+(
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::VirtualExpression<cla3p::csc::XxMatrix<T_Int,T_Scalar>, T_Right>& B)
{
	return (A.virtualize() + B);
}

/*-------------------------------------------------*/

//
// Virtual + XxObject
//

template <typename T_Scalar, typename T_Left>
cla3p::VirtualPlus<
	cla3p::dns::XxVector<T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::dns::XxVector<T_Scalar>>>
operator+(
		const cla3p::VirtualExpression<cla3p::dns::XxVector<T_Scalar>, T_Left>& A,
		const cla3p::dns::XxVector<T_Scalar>& B)
{
	return (A + B.virtualize());
}

template <typename T_Scalar, typename T_Left>
cla3p::VirtualPlus<
	cla3p::dns::XxMatrix<T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>
operator+(
		const cla3p::VirtualExpression<cla3p::dns::XxMatrix<T_Scalar>, T_Left>& A,
		const cla3p::dns::XxMatrix<T_Scalar>& B)
{
	return (A + B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
cla3p::VirtualPlus<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>
operator+(
		const cla3p::VirtualExpression<cla3p::csc::XxMatrix<T_Int,T_Scalar>, T_Left>& A,
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B)
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
cla3p::VirtualMinus<T_Result,T_Left,T_Right>
operator-(
	const cla3p::VirtualExpression<T_Result, T_Left>& A, 
	const cla3p::VirtualExpression<T_Result, T_Right>& B)
{
	return cla3p::VirtualMinus<T_Result,T_Left,T_Right>(A, B);
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
cla3p::alias::VirtualMinus_vec<T_Scalar>
operator-(
		const cla3p::dns::XxVector<T_Scalar>& x,
		const cla3p::dns::XxVector<T_Scalar>& y)
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
cla3p::alias::VirtualMinus_dns<T_Scalar>
operator-(
		const cla3p::dns::XxMatrix<T_Scalar>& A,
		const cla3p::dns::XxMatrix<T_Scalar>& B)
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
cla3p::alias::VirtualMinus_csc<T_Int,T_Scalar>
operator-(
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A.virtualize() - B.virtualize());
}

/*-------------------------------------------------*/

//
// XxObject - Virtual
//

template <typename T_Scalar, typename T_Right>
cla3p::VirtualMinus<
	cla3p::dns::XxVector<T_Scalar>,
	cla3p::VirtualObject<cla3p::dns::XxVector<T_Scalar>>,
	T_Right>
operator-(
		const cla3p::dns::XxVector<T_Scalar>& A,
		const cla3p::VirtualExpression<cla3p::dns::XxVector<T_Scalar>, T_Right>& B)
{
	return (A.virtualize() - B);
}

template <typename T_Scalar, typename T_Right>
cla3p::VirtualMinus<
	cla3p::dns::XxMatrix<T_Scalar>,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>,
	T_Right>
operator-(
		const cla3p::dns::XxMatrix<T_Scalar>& A,
		const cla3p::VirtualExpression<cla3p::dns::XxMatrix<T_Scalar>, T_Right>& B)
{
	return (A.virtualize() - B);
}

template <typename T_Int, typename T_Scalar, typename T_Right>
cla3p::VirtualMinus<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>,
	T_Right>
operator-(
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A,
		const cla3p::VirtualExpression<cla3p::csc::XxMatrix<T_Int,T_Scalar>, T_Right>& B)
{
	return (A.virtualize() - B);
}

/*-------------------------------------------------*/

//
// Virtual + XxObject
//

template <typename T_Scalar, typename T_Left>
cla3p::VirtualMinus<
	cla3p::dns::XxVector<T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::dns::XxVector<T_Scalar>>>
operator-(
		const cla3p::VirtualExpression<cla3p::dns::XxVector<T_Scalar>, T_Left>& A,
		const cla3p::dns::XxVector<T_Scalar>& B)
{
	return (A - B.virtualize());
}

template <typename T_Scalar, typename T_Left>
cla3p::VirtualMinus<
	cla3p::dns::XxMatrix<T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>
operator-(
		const cla3p::VirtualExpression<cla3p::dns::XxMatrix<T_Scalar>, T_Left>& A,
		const cla3p::dns::XxMatrix<T_Scalar>& B)
{
	return (A - B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
cla3p::VirtualMinus<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>
operator-(
		const cla3p::VirtualExpression<cla3p::csc::XxMatrix<T_Int,T_Scalar>, T_Left>& A,
		const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B)
{
	return (A - B.virtualize());
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_ADD_HPP_
