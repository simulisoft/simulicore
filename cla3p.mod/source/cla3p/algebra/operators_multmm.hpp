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

#ifndef CLA3P_OPERATORS_MULTMM_HPP_
#define CLA3P_OPERATORS_MULTMM_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/virtuals/virtual_product.hpp"

/*-------------------------------------------------*/

/*
 * Virtual x XxMatrix
 */

template <typename T_Scalar, typename T_Left>
cla3p::VirtualProduct<
	cla3p::dns::XxMatrix<T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>> 
operator*(
	const cla3p::VirtualExpression<typename T_Left::result_type, T_Left>& left, 
	const cla3p::dns::XxMatrix<T_Scalar>& B) 
{ 
	return cla3p::VirtualProduct<
		cla3p::dns::XxMatrix<T_Scalar>,
		T_Left,
		cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>(left.self(), B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
cla3p::VirtualProduct<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	T_Left, 
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>
operator*(
	const cla3p::VirtualExpression<typename T_Left::result_type, T_Left>& left, 
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B) 
{ 
	return cla3p::VirtualProduct<
		cla3p::csc::XxMatrix<T_Int,T_Scalar>,
		T_Left,
		cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>(left.self(), B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Multiplies a matrix with a matrix.
 * @details Performs the operation <b>A * B</b>
 * @param[in] A The lhs input matrix.
 * @param[in] B The rhs input matrix.
 * @return The virtual product.
 */
template <typename T_Scalar>
cla3p::VirtualProduct<
	cla3p::dns::XxMatrix<T_Scalar>,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>> 
operator*(
	const cla3p::dns::XxMatrix<T_Scalar>& A, 
	const cla3p::dns::XxMatrix<T_Scalar>& B) 
{
	return cla3p::VirtualProduct<
		cla3p::dns::XxMatrix<T_Scalar>,
		cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>,
		cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>(A.virtualize(), B.virtualize()); 
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Multiplies a sparse matrix with a dense matrix.
 * @details Performs the operation <b>A * B</b>
 * @warning The virtual matrix product `A * B.transpose()` 
 *          where A is sparse and B is dense calculates the transpose of B explicitly.@n
 *          If you plan to use the specific product in a loop, 
 *          it is recommended to pre-calculate the transpose of B.
 * @param[in] A The lhs input matrix.
 * @param[in] B The rhs input matrix.
 * @return The resulting dense matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::VirtualProduct<
	cla3p::dns::XxMatrix<T_Scalar>,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>, 
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>
operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::dns::XxMatrix<T_Scalar>& B) 
{
	return cla3p::VirtualProduct<
		cla3p::dns::XxMatrix<T_Scalar>,
		cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>,
		cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>(A.virtualize(), B.virtualize()); 
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Multiplies a sparse matrix with a sparse matrix.
 * @details Performs the operation <b>A * B</b>
 * @param[in] A The lhs input matrix.
 * @param[in] B The rhs input matrix.
 * @return The resulting sparse matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::VirtualProduct<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>> 
operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B) 
{
	return cla3p::VirtualProduct<
		cla3p::csc::XxMatrix<T_Int,T_Scalar>,
		cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>,
		cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>(A.virtualize(), B.virtualize()); 
}

#endif // CLA3P_OPERATORS_MULTMM_HPP_
