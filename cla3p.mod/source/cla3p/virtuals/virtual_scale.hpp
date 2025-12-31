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

#ifndef CLA3P_VIRTUAL_SCALE_HPP_
#define CLA3P_VIRTUAL_SCALE_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/virtuals/virtual_object.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The virtual scaling expression class.
 * @details Virtually scales expressions.
 */
template <typename T_Result, typename T_Virtual>
class VirtualScale : public VirtualExpression<T_Result, VirtualScale<T_Result, T_Virtual>> {

	private:
		using T_Scalar = typename T_Result::value_type;

	public:
		explicit VirtualScale(const VirtualExpression<T_Result, T_Virtual>& src, const T_Scalar& coeff)
		: m_src(src.self()), m_coeff(coeff) { }

		~VirtualScale() {}

		void evaluateOnNew(T_Result& dest) const override
		{
			dest.clear();
			m_src.evaluateOnNew(dest);
			dest.iscale(m_coeff);
		}

		void evaluateOnExisting(T_Result& dest) const override
		{
			m_src.evaluateOnExisting(dest);
			dest.iscale(m_coeff);
		}

		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override
		{
			m_src.accumulateOnExisting(dest, coeff * m_coeff);
		}

		const T_Virtual& get() const { return m_src; }
		const T_Scalar& coeff() const { return m_coeff; }

	private:
		T_Virtual m_src;
		T_Scalar m_coeff;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_SCALE_HPP_
