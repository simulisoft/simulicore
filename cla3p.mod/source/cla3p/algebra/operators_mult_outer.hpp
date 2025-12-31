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

#ifndef CLA3P_OPERATORS_MULT_OUTER_HPP_
#define CLA3P_OPERATORS_MULT_OUTER_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/virtuals/virtual_rowvec.hpp"
#include "cla3p/virtuals/virtual_outer.hpp"

/*-------------------------------------------------*/

/*
 * XxVector x VirtualRowvec
 */
template <typename T_Scalar>
cla3p::VirtualOuter<
	cla3p::dns::XxMatrix<T_Scalar>,
	cla3p::VirtualObject<cla3p::dns::XxVector<T_Scalar>>>
operator*(
		const cla3p::dns::XxVector<T_Scalar>& X,
		const cla3p::VirtualRowvec<T_Scalar>& vY)
{
	return cla3p::VirtualOuter<
		cla3p::dns::XxMatrix<T_Scalar>,
		cla3p::VirtualObject<cla3p::dns::XxVector<T_Scalar>>>(X.virtualize(), vY);
}

/*
 * VirtualExpression x VirtualRowvec
 */
template <typename T_Scalar, typename T_Virtual>

cla3p::VirtualOuter<cla3p::dns::XxMatrix<T_Scalar>,T_Virtual>
operator*(
		const cla3p::VirtualExpression<cla3p::dns::XxVector<T_Scalar>, T_Virtual>& vX,
		const cla3p::VirtualRowvec<T_Scalar>& vY)
{
	return cla3p::VirtualOuter<cla3p::dns::XxMatrix<T_Scalar>,T_Virtual>(vX, vY);
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_MULT_OUTER_HPP_
