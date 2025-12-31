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

#ifndef CLA3P_OPERATORS_MULT_INNER_HPP_
#define CLA3P_OPERATORS_MULT_INNER_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/virtuals/virtual_rowvec.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace dns { template <typename T_Scalar> class XxVector; }
} // namespace cla3p
/*-------------------------------------------------*/

/*
 * VirtualRowvec x XxVector
 */
template <typename T_Scalar>
T_Scalar operator*(
		const cla3p::VirtualRowvec<T_Scalar>& vX,
		const cla3p::dns::XxVector<T_Scalar>& Y)
{
	return vX.evaluateInner(Y);
}

/*
 * VirtualRowvec x VirtualExpression
 */
template <typename T_Vector, typename T_Virtual>
typename T_Vector::value_type operator*(
		const cla3p::VirtualRowvec<typename T_Vector::value_type>& vX,
		const cla3p::VirtualExpression<T_Vector, T_Virtual>& vY)
{
	return (vX * vY.evaluate());
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_MULT_INNER_HPP_
