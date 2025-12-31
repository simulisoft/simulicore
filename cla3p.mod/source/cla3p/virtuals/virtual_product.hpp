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

#ifndef CLA3P_VIRTUAL_PRODUCT_HPP_
#define CLA3P_VIRTUAL_PRODUCT_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The virtual product class.
 * @details Virtually multiplies two virtual expressions.
 */
template <typename T_Result, typename T_Left, typename T_Right>
class VirtualProduct : public VirtualExpression<T_Result, VirtualProduct<T_Result, T_Left, T_Right>> {

	private:
		using T_Scalar = typename T_Result::value_type;

	public:
		VirtualProduct(
			const VirtualExpression<typename T_Left::result_type,T_Left>& left, 
			const VirtualExpression<typename T_Right::result_type,T_Right>& right)
		: m_left(left.self()), m_right(right.self()) {}
		~VirtualProduct() {}

		void evaluateOnNew(T_Result& dest) const override;
		void evaluateOnExisting(T_Result& dest) const override;
		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override;

	private:
		T_Left m_left;
		T_Right m_right;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#include "cla3p/virtuals/virtual_product_specs_impl.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualProduct<T_Result, T_Left, T_Right>::evaluateOnNew(T_Result& dest) const
{ 
	VirtualProductEvaluateOnNewSpec(m_left, m_right, dest); 
}
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualProduct<T_Result, T_Left, T_Right>::evaluateOnExisting(T_Result& dest) const
{
	VirtualProductEvaluateOnExistingSpec(m_left, m_right, dest);
}
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualProduct<T_Result, T_Left, T_Right>::accumulateOnExisting(T_Result& dest, T_Scalar coeff) const
{
	VirtualProductAccumulateOnExistingSpec(m_left, m_right, dest, coeff);
}
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_PRODUCT_HPP_
