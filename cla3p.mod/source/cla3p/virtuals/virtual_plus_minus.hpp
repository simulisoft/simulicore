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

#ifndef CLA3P_VIRTUAL_PLUS_MINUS_HPP_
#define CLA3P_VIRTUAL_PLUS_MINUS_HPP_

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
 * @brief The virtual addition expression class.
 * @details Virtually adds two virtual expressions.
 */
template <typename T_Result, typename T_Left, typename T_Right>
class VirtualPlus : public VirtualExpression<T_Result, VirtualPlus<T_Result, T_Left, T_Right>> {

	private:
		using T_Scalar = typename T_Result::value_type;

	public:
		VirtualPlus(const VirtualExpression<T_Result,T_Left>& left, const VirtualExpression<T_Result,T_Right>& right) 
		: m_left(left.self()), m_right(right.self()) {}
		~VirtualPlus() {}

		void evaluateOnNew(T_Result& dest) const override
		{
			dest.clear();
			m_left.evaluateOnNew(dest);
			m_right.accumulateOnExisting(dest, T_Scalar(1));
		}

		void evaluateOnExisting(T_Result& dest) const override;

		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override
		{
			m_left.accumulateOnExisting(dest, coeff);
			m_right.accumulateOnExisting(dest, coeff);
		}

	private:
		T_Left m_left;
		T_Right m_right;
};

/*-------------------------------------------------*/
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualPlusEvaluateOnExistingSpec(
	const VirtualExpression<typename T_Left::result_type,T_Left>& left, 
	const VirtualExpression<typename T_Right::result_type,T_Right>& right, 
	T_Result& dest)
{ 
	left.evaluateOnExisting(dest);
	right.accumulateOnExisting(dest, typename T_Result::value_type(1));
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Left, typename T_Right>
void VirtualPlusEvaluateOnExistingSpec(
	const VirtualExpression<typename T_Left::result_type,T_Left>& left, 
	const VirtualExpression<typename T_Right::result_type,T_Right>& right, 
	csc::XxMatrix<T_Int,typename T_Right::return_type::value_type>& dest)
{ 
	csc::XxMatrix<T_Int,typename T_Right::return_type::value_type> tmp;
	evaluateOnNew(tmp);
	dest = tmp;
}
/*-------------------------------------------------*/
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualPlus<T_Result, T_Left, T_Right>::evaluateOnExisting(T_Result& dest) const
{
	VirtualPlusEvaluateOnExistingSpec<T_Result, T_Left, T_Right>(m_left, m_right, dest); 
}
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The virtual subtraction expression class.
 * @details Virtually subtracts two expressions.
 */
template <typename T_Result, typename T_Left, typename T_Right>
class VirtualMinus : public VirtualExpression<T_Result, VirtualMinus<T_Result, T_Left, T_Right>> {

	private:
		using T_Scalar = typename T_Result::value_type;

	public:
		VirtualMinus(const VirtualExpression<T_Result,T_Left>& left, const VirtualExpression<T_Result,T_Right>& right) 
		: m_left(left.self()), m_right(right.self()) {}
		~VirtualMinus() {}

		void evaluateOnNew(T_Result& dest) const override
		{
			dest.clear();
			m_left.evaluateOnNew(dest);
			m_right.accumulateOnExisting(dest, T_Scalar(-1));
		}

		void evaluateOnExisting(T_Result& dest) const override;

		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override
		{
			m_left.accumulateOnExisting(dest, coeff);
			m_right.accumulateOnExisting(dest, -coeff); 
		}

	private:
		T_Left m_left;
		T_Right m_right;
};

/*-------------------------------------------------*/
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualMinusEvaluateOnExistingSpec(
	const VirtualExpression<typename T_Left::result_type,T_Left>& left, 
	const VirtualExpression<typename T_Right::result_type,T_Right>& right, 
	T_Result& dest)
{ 
	left.evaluateOnExisting(dest);
	right.accumulateOnExisting(dest, typename T_Result::value_type(-1));
}
/*-------------------------------------------------*/
template <typename T_Int, typename T_Left, typename T_Right>
void VirtualMinusEvaluateOnExistingSpec(
	const VirtualExpression<typename T_Left::result_type,T_Left>& left, 
	const VirtualExpression<typename T_Right::result_type,T_Right>& right, 
	csc::XxMatrix<T_Int,typename T_Right::return_type::value_type>& dest)
{ 
	csc::XxMatrix<T_Int,typename T_Right::return_type::value_type> tmp;
	evaluateOnNew(tmp);
	dest = tmp;
}
/*-------------------------------------------------*/
template <typename T_Result, typename T_Left, typename T_Right>
void VirtualMinus<T_Result, T_Left, T_Right>::evaluateOnExisting(T_Result& dest) const
{
	VirtualMinusEvaluateOnExistingSpec<T_Result, T_Left, T_Right>(m_left, m_right, dest); 
}
/*-------------------------------------------------*/

namespace alias { 

template <typename T_Scalar>
using VirtualPlus_vec = VirtualPlus<
	dns::XxVector<T_Scalar>,
	VirtualObj_vec<T_Scalar>,
	VirtualObj_vec<T_Scalar>>;

template <typename T_Scalar>
using VirtualPlus_dns = VirtualPlus<
	dns::XxMatrix<T_Scalar>,
	VirtualObj_dns<T_Scalar>,
	VirtualObj_dns<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualPlus_csc = VirtualPlus<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>>;

template <typename T_Scalar>
using VirtualMinus_vec = VirtualMinus<
	dns::XxVector<T_Scalar>,
	VirtualObj_vec<T_Scalar>,
	VirtualObj_vec<T_Scalar>>;

template <typename T_Scalar>
using VirtualMinus_dns = VirtualMinus<
	dns::XxMatrix<T_Scalar>,
	VirtualObj_dns<T_Scalar>,
	VirtualObj_dns<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualMinus_csc = VirtualMinus<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>>;

} // namespace alias

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_PLUS_MINUS_HPP_
