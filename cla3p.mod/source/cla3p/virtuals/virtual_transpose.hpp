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

#ifndef CLA3P_VIRTUAL_TRANSPOSE_HPP_
#define CLA3P_VIRTUAL_TRANSPOSE_HPP_

/**
 * @file
 */
#include "cla3p/types/enums.hpp"
#include "cla3p/virtuals/virtual_expression.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }

/*-------------------------------------------------*/

template <typename T_Scalar>
void VirtualTransposeEvaluateOnNewSpec(const dns::XxMatrix<T_Scalar>& src, bool conj, dns::XxMatrix<T_Scalar>& dest);
template <typename T_Scalar>
void VirtualTransposeEvaluateOnExistingSpec(const dns::XxMatrix<T_Scalar>& src, bool conj, dns::XxMatrix<T_Scalar>& dest);
template <typename T_Scalar>
void VirtualTransposeAccumulateOnExistingSpec(const dns::XxMatrix<T_Scalar>& src, bool conj, dns::XxMatrix<T_Scalar>& dest, T_Scalar coeff);

template <typename T_Int, typename T_Scalar>
void VirtualTransposeEvaluateOnNewSpec(const csc::XxMatrix<T_Int,T_Scalar>& src, bool conj, csc::XxMatrix<T_Int, T_Scalar>& dest);
template <typename T_Int, typename T_Scalar>
void VirtualTransposeEvaluateOnExistingSpec(const csc::XxMatrix<T_Int,T_Scalar>& src, bool conj, csc::XxMatrix<T_Int, T_Scalar>& dest);
template <typename T_Int, typename T_Scalar>
void VirtualTransposeAccumulateOnExistingSpec(const csc::XxMatrix<T_Int,T_Scalar>& src, bool conj, csc::XxMatrix<T_Int, T_Scalar>& dest, T_Scalar coeff);

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The virtual (conjugate-)transpose class.
 */
template <typename T_Result>
class VirtualTranspose : public VirtualExpression<T_Result, VirtualTranspose<T_Result>> {

	private:
		using T_Scalar = typename T_Result::value_type;

	public:
		explicit VirtualTranspose(const T_Result& src, bool conj) : m_src(src), m_conj(conj) {}
		~VirtualTranspose() {}

		void evaluateOnNew(T_Result& dest) const override { VirtualTransposeEvaluateOnNewSpec(m_src, m_conj, dest); }
		void evaluateOnExisting(T_Result& dest) const override { VirtualTransposeEvaluateOnExistingSpec(m_src, m_conj, dest); }
		void accumulateOnExisting(T_Result& dest, T_Scalar coeff) const override { VirtualTransposeAccumulateOnExistingSpec(m_src, m_conj, dest, coeff); }

		const T_Result& get() const { return m_src; }
		op_t op() const { return m_conj ? op_t::C : op_t::T; }

	private:
		const T_Result& m_src;
		bool m_conj;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_TRANSPOSE_HPP_
