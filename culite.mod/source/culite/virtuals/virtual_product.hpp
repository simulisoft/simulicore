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

#ifndef CULITE_VIRTUAL_PRODUCT_HPP_
#define CULITE_VIRTUAL_PRODUCT_HPP_

/**
 * @file
 */

#include "culite/virtuals/virtual_expression.hpp"
#include "culite/virtuals/virtual_object.hpp"

/*-------------------------------------------------*/
namespace culite { 
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

namespace alias { 

template <typename T_Scalar>
using VirtualProd_dnsmv = VirtualProduct<
	dns::XxVector<T_Scalar>,
	VirtualObj_dns<T_Scalar>,
	VirtualObj_vec<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualProd_csrmv = VirtualProduct<
	dns::XxVector<T_Scalar>,
	VirtualObj_csr<T_Int,T_Scalar>,
	VirtualObj_vec<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualProd_cscmv = VirtualProduct<
	dns::XxVector<T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>,
	VirtualObj_vec<T_Scalar>>;

template <typename T_Scalar>
using VirtualProd_dnsmm = VirtualProduct<
	dns::XxMatrix<T_Scalar>,
	VirtualObj_dns<T_Scalar>,
	VirtualObj_dns<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualProd_csrmm = VirtualProduct<
	dns::XxMatrix<T_Scalar>,
	VirtualObj_csr<T_Int,T_Scalar>,
	VirtualObj_dns<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualProd_cscmm = VirtualProduct<
	dns::XxMatrix<T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>,
	VirtualObj_dns<T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualProd_csrspmm = VirtualProduct<
	csr::XxMatrix<T_Int,T_Scalar>,
	VirtualObj_csr<T_Int,T_Scalar>,
	VirtualObj_csr<T_Int,T_Scalar>>;

template <typename T_Int, typename T_Scalar>
using VirtualProd_cscspmm = VirtualProduct<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>,
	VirtualObj_csc<T_Int,T_Scalar>>;

} // namespace alias

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#include "culite/virtuals/virtual_product_specs_impl.hpp"

/*-------------------------------------------------*/
namespace culite { 
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
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_VIRTUAL_PRODUCT_HPP_
