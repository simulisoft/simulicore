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

#ifndef CLA3P_DNS_XXVECTOR_HPP_
#define CLA3P_DNS_XXVECTOR_HPP_

/**
 * @file
 */

#include <string>
#include <ostream>

#include "cla3p/generic/guard.hpp"
#include "cla3p/dense/dns_xivector.hpp"

#include "cla3p/virtuals/virtual_rowvec.hpp"
#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_transpose.hpp"
#include "cla3p/virtuals/virtual_conjugate.hpp"
#include "cla3p/virtuals/virtual_scale.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace prm { template <typename T_Int> class PxMatrix; }

/*-------------------------------------------------*/
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The floating point dense vector class.
 */
template <typename T_Scalar>
class XxVector : public XiVector<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		explicit XxVector(const XiVector<T_Scalar>& other);
		XxVector(XiVector<T_Scalar>&& other);

		XxVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);
		XxVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

		template <typename T_Virtual>
		explicit XxVector(const VirtualExpression<XxVector<T_Scalar>,T_Virtual>& v) { operator=(v); }
		template <typename T_Virtual>
		XxVector(VirtualExpression<XxVector<T_Scalar>,T_Virtual>&& v) { operator=(std::move(v)); }

		template <typename T_Virtual>
		XxVector<T_Scalar>& operator=(const VirtualExpression<XxVector<T_Scalar>,T_Virtual>& v) { evaluateFrom(v); return *this; }
		template <typename T_Virtual>
		XxVector<T_Scalar>& operator=(VirtualExpression<XxVector<T_Scalar>,T_Virtual>&& v) { evaluateFrom(v); return *this; }

		VirtualObject<XxVector<T_Scalar>> virtualize() const { return VirtualObject<XxVector<T_Scalar>>(*this); }

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_vector_docs::constructor()
		 */
		XxVector();

		/**
		 * @copydoc standard_vector_docs::dim_constructor()
		 */
		explicit XxVector(int_t n);

		/**
		 * @copydoc standard_vector_docs::aux_constructor()
		 */
		explicit XxVector(int_t n, T_Scalar *vals, bool bind);

		/**
		 * @copydoc standard_docs::copy_constructor()
		 */
		XxVector(const XxVector<T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_constructor()
		 */
		XxVector(XxVector<T_Scalar>&& other) = default;

		/**
		 * @copydoc standard_vector_docs::destructor()
		 */
		~XxVector();

		/** @} */

		/** 
		 * @name Operators
		 * @{
		 */

		/**
		 * @copydoc standard_docs::copy_assignment()
		 */
		XxVector<T_Scalar>& operator=(const XxVector<T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_assignment()
		 */
		XxVector<T_Scalar>& operator=(XxVector<T_Scalar>&& other) = default;

		/**
		 * @copydoc standard_matrix_docs::fill()
		 */
		void operator=(T_Scalar val);

		/**
		 * @copydoc standard_docs::virtual_negate_operator()
		 */
		VirtualScale<XxVector<T_Scalar>,VirtualObject<XxVector<T_Scalar>>> operator-() const;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @copydoc standard_docs::iscale()
		 */
		void iscale(T_Scalar val);

		/**
		 * @copydoc standard_vector_docs::virtual_transpose()
		 */
		VirtualRowvec<T_Scalar> transpose() const;

		/**
		 * @copydoc standard_vector_docs::virtual_ctranspose()
		 */
		VirtualRowvec<T_Scalar> ctranspose() const;

		/**
		 * @copydoc standard_vector_docs::virtual_conjugate()
		 */
		VirtualConjugate<XxVector<T_Scalar>> conjugate() const;

		/**
		 * @copydoc standard_vector_docs::iconjugate()
		 */
		void iconjugate();

		/**
		 * @copydoc standard_docs::normOne()
		 */
		T_RScalar normOne() const;

		/**
		 * @copydoc standard_docs::normInf()
		 */
		T_RScalar normInf() const;

		/**
		 * @copydoc standard_docs::normEuc()
		 */
		T_RScalar normEuc() const;

		/**
		 * @copydoc standard_vector_docs::permute_left()
		 */
		XxVector<T_Scalar> permuteLeft(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @copydoc standard_vector_docs::permute_left_dst()
		 */
		void permuteLeft(const prm::PxMatrix<int_t>& P, XxVector<T_Scalar>& dest) const;

		/**
		 * @copydoc standard_vector_docs::rblock()
		 */
		XxVector<T_Scalar> block(int_t ibgn, int_t ni) const;

		/**
		 * @copydoc standard_vector_docs::rblock()
		 */
		XxVector<T_Scalar> rblock(int_t ibgn, int_t ni);

		/**
		 * @copydoc standard_vector_docs::rblock()
		 */
		Guard<XxVector<T_Scalar>> rblock(int_t ibgn, int_t ni) const;

		/**
		 * @copydoc standard_vector_docs::set_block()
		 */
		void setBlock(int_t ibgn, const XxVector<T_Scalar>& src);

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @copydoc standard_vector_docs::random()
		 */
		static XxVector<T_Scalar> random(int_t n, T_RScalar lo = T_RScalar(0), T_RScalar hi = T_RScalar(1));

		/** @} */

	protected:
		template <typename T_Virtual>
		void evaluateFrom(const VirtualExpression<XxVector<T_Scalar>,T_Virtual>& v)
		{
			if(*this) {
				v.evaluateOnExisting(*this);
			} else {
				v.evaluateOnNew(*this);
			}
		}

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DNS_XXVECTOR_HPP_
