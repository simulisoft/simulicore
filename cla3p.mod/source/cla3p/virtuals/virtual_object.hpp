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

#ifndef CLA3P_VIRTUAL_OBJECT_HPP_
#define CLA3P_VIRTUAL_OBJECT_HPP_

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
 * @brief The virtual object class.
 * @details Virtualizes a non-virtual object, making it compatible with all virtual expressions.
 */
template <typename T_Result>
class VirtualObject : public VirtualExpression<T_Result, VirtualObject<T_Result>> {

	private:
		using T_Scalar = typename T_Result::value_type;
		
	public:
		explicit VirtualObject(const T_Result& obj) : m_obj(obj) {}
		~VirtualObject() {}

		void evaluateOnNew(T_Result& dest) const override
		{ 
			dest.clear(); 
			dest = m_obj; 
		}

		void evaluateOnExisting(T_Result& dest) const override
		{ 
			dest = m_obj; 
		}

		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override
		{ 
			ops::update(coeff, m_obj, dest); 
		};

		const T_Result& get() const { return m_obj; }

	private:
		const T_Result& m_obj;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_OBJECT_HPP_