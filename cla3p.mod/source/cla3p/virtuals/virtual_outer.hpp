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

#ifndef CLA3P_VIRTUAL_OUTER_HPP_
#define CLA3P_VIRTUAL_OUTER_HPP_

/**
 * @file
 */

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_rowvec.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }

/**
 * @nosubgrouping
 * @brief The virtual outer product expression class.
 * @details Virtually computes outer products of expression with VirtualRowvec.
 */
template <typename T_Result, typename T_Virtual>
class VirtualOuter : public VirtualExpression<T_Result, VirtualOuter<T_Result, T_Virtual>> {

	private:
		using T_Scalar = typename T_Result::value_type;
		using T_Vector = dns::XxVector<T_Scalar>;

	public:
		explicit VirtualOuter(const VirtualExpression<T_Vector, T_Virtual>& left, const VirtualRowvec<T_Scalar>& right)
		: m_left(left.self()), m_right(right) {}

		~VirtualOuter() {}

		void evaluateOnNew(T_Result& dest) const override;
		void evaluateOnExisting(T_Result& dest) const override;
		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override;

	private:
		T_Virtual m_left;
		const VirtualRowvec<T_Scalar>& m_right;
};

/*-------------------------------------------------*/

template <typename T_Scalar>
void VirtualOuterEvaluateOnNewSpec(
	const VirtualObject<dns::XxVector<T_Scalar>>& left, 
	const VirtualRowvec<T_Scalar>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	right.evaluateOuterOnNew(T_Scalar(1), left.get(), dest);
}
/*-------------------------------------------------*/
template <typename T_Virtual>
void VirtualOuterEvaluateOnNewSpec(
	const VirtualExpression<typename T_Virtual::result_type, T_Virtual>& left, 
	const VirtualRowvec<typename T_Virtual::result_type::value_type>& right, 
	dns::XxMatrix<typename T_Virtual::result_type::value_type>& dest)
{
	right.evaluateOuterOnNew(typename T_Virtual::result_type::value_type(1), left.evaluate(), dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualOuterEvaluateOnExistingSpec(
	const VirtualObject<dns::XxVector<T_Scalar>>& left, 
	const VirtualRowvec<T_Scalar>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	right.evaluateOuterOnExisting(T_Scalar(1), left.get(), dest);
}
/*-------------------------------------------------*/
template <typename T_Virtual>
void VirtualOuterEvaluateOnExistingSpec(
	const VirtualExpression<typename T_Virtual::result_type, T_Virtual>& left, 
	const VirtualRowvec<typename T_Virtual::result_type::value_type>& right, 
	dns::XxMatrix<typename T_Virtual::result_type::value_type>& dest)
{
	right.evaluateOuterOnExisting(typename T_Virtual::result_type::value_type(1), left.evaluate(), dest);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void VirtualOuterAccumulateOnExistingSpec(
	const VirtualObject<dns::XxVector<T_Scalar>>& left, 
	const VirtualRowvec<T_Scalar>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	right.accumulateOuterOnExisting(coeff, left.get(), dest);
}
/*-------------------------------------------------*/
template <typename T_Virtual>
void VirtualOuterAccumulateOnExistingSpec(
	const VirtualExpression<typename T_Virtual::result_type, T_Virtual>& left, 
	const VirtualRowvec<typename T_Virtual::result_type::value_type>& right, 
	dns::XxMatrix<typename T_Virtual::result_type::value_type>& dest, 
	typename T_Virtual::result_type::value_type coeff)
{
	right.accumulateOuterOnExisting(coeff, left.evaluate(), dest);
}
/*-------------------------------------------------*/
template <typename T_Result, typename T_Virtual>
void VirtualOuter<T_Result, T_Virtual>::evaluateOnNew(T_Result& dest) const
{
	VirtualOuterEvaluateOnNewSpec(m_left, m_right, dest);
}
/*-------------------------------------------------*/
template <typename T_Result, typename T_Virtual>
void VirtualOuter<T_Result, T_Virtual>::evaluateOnExisting(T_Result& dest) const
{
	VirtualOuterEvaluateOnExistingSpec(m_left, m_right, dest);
}
/*-------------------------------------------------*/
template <typename T_Result, typename T_Virtual>
void VirtualOuter<T_Result, T_Virtual>::accumulateOnExisting(T_Result& dest, T_Scalar coeff) const
{
	VirtualOuterAccumulateOnExistingSpec(m_left, m_right, dest, coeff);
}
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_OUTER_HPP_
