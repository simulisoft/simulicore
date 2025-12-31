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

#ifndef CLA3P_VIRTUAL_CONJUGATE_HPP_
#define CLA3P_VIRTUAL_CONJUGATE_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/algebra/functional_update.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The virtual conjugate class.
 */
template <typename T_Result>
class VirtualConjugate : public VirtualExpression<T_Result, VirtualConjugate<T_Result>> {

	private:
		using T_Scalar = typename T_Result::value_type;

	public:
		explicit VirtualConjugate(const T_Result& src) : m_src(src) {}
		~VirtualConjugate() {}

		void evaluateOnNew(T_Result& dest) const override 
		{ 
			dest.clear();
			dest = m_src;
			dest.iconjugate();
		}

		void evaluateOnExisting(T_Result& dest) const override
		{
			dest = m_src;
			dest.iconjugate();
		}

		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override
		{ 
			T_Result tmp;
			evaluateOnNew(tmp);
			ops::update(coeff, tmp, dest);
		}


	private:
		const T_Result& m_src;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_CONJUGATE_HPP_
