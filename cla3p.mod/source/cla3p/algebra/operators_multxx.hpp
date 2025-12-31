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

#ifndef CLA3P_OPERATORS_MULTXX_HPP_
#define CLA3P_OPERATORS_MULTXX_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_product.hpp"

/*-------------------------------------------------*/

/*
 * Virtual x Virtual
 */

template <typename T_Lhs, typename T_Rhs>
cla3p::VirtualProduct<typename T_Rhs::result_type, T_Lhs, T_Rhs>
operator*(
	const cla3p::VirtualExpression<typename T_Lhs::result_type, T_Lhs>& A, 
	const cla3p::VirtualExpression<typename T_Rhs::result_type, T_Rhs>& B)
{
	return cla3p::VirtualProduct<typename T_Rhs::result_type,T_Lhs,T_Rhs>(A.self(), B.self());
}

/*
 * XxMatrix x Virtual
 */

template <typename T_Scalar, typename T_Right>
cla3p::VirtualProduct<
	typename T_Right::result_type,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>, 
	T_Right>
operator*(
	const cla3p::dns::XxMatrix<T_Scalar>& A, 
	const cla3p::VirtualExpression<typename T_Right::result_type, T_Right>& right) 
{ 
	return cla3p::VirtualProduct<
		typename T_Right::result_type,
		cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>,
		T_Right>(A.virtualize(), right.self());
}

template <typename T_Int, typename T_Scalar, typename T_Right>
cla3p::VirtualProduct<
	typename T_Right::result_type,
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>, 
	T_Right>
operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::VirtualExpression<typename T_Right::result_type, T_Right>& right)
{ 
	return cla3p::VirtualProduct<
		typename T_Right::result_type,
		cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>,
		T_Right>(A.virtualize(), right.self());
}

#endif // CLA3P_OPERATORS_MULTXX_HPP_
