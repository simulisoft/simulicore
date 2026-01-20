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

#ifndef CLA3P_OPERATORS_MULTMV_HPP_
#define CLA3P_OPERATORS_MULTMV_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

#include "cla3p/virtuals/virtual_product.hpp"

/*-------------------------------------------------*/

/*
 * Virtual x XxVector
 */
template <typename T_Left>
cla3p::VirtualProduct<
	cla3p::dns::XxVector<typename T_Left::result_type::value_type>,
	T_Left,
	cla3p::alias::VirtualObj_vec<typename T_Left::result_type::value_type>> 
operator*(
	const cla3p::VirtualExpression<typename T_Left::result_type, T_Left>& left, 
	const cla3p::dns::XxVector<typename T_Left::result_type::value_type>& x) 
{ 
	return cla3p::VirtualProduct<
		cla3p::dns::XxVector<typename T_Left::result_type::value_type>,
		T_Left,
		cla3p::alias::VirtualObj_vec<typename T_Left::result_type::value_type>>(left.self(), x.virtualize());
}

/*
 * XxMatrix x XxVector
 */

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Matrix-vector multiplication for dense matrices.
 * @details Performs the operation @f$ A \cdot x @f$.
 * @param[in] A The input dense matrix.
 * @param[in] x The input vector.
 * @return The virtual product representing the multiplication.
 */
template <typename T_Scalar>
cla3p::alias::VirtualProd_dnsmv<T_Scalar>
operator*(
	const cla3p::dns::XxMatrix<T_Scalar>& A, 
	const cla3p::dns::XxVector<T_Scalar>& x) 
{ 
	return cla3p::alias::VirtualProd_dnsmv<T_Scalar>(A.virtualize(), x.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Matrix-vector multiplication for sparse matrices.
 * @details Performs the operation @f$ A \cdot x @f$.
 * @param[in] A The input sparse matrix.
 * @param[in] x The input vector.
 * @return The virtual product representing the multiplication.
 */
template <typename T_Int, typename T_Scalar>
cla3p::alias::VirtualProd_cscmv<T_Int,T_Scalar>
operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::dns::XxVector<T_Scalar>& x) 
{ 
	return cla3p::alias::VirtualProd_cscmv<T_Int,T_Scalar>(A.virtualize(), x.virtualize());
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_MULTMV_HPP_
