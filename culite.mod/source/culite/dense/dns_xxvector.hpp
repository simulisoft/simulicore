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

#ifndef CULITE_DNS_XXVECTOR_HPP_
#define CULITE_DNS_XXVECTOR_HPP_

/**
 * @file
 */

#include <string>
#include <ostream>

#include "culite/dense/dns_xivector.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The floating point device dense vector class.
 * @details Represents a floating-point dense vector stored on the GPU device.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
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
		 * @copydoc standard_docs::virtual_negate_operator()
		 */
		XxVector<T_Scalar> operator-() const; // TODO: use virtuals

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @copydoc standard_docs::iscale()
		 */
		void iscale(const T_Scalar& val);

		/* TODO: use virtuals
		 * @copydoc standard_vector_docs::virtual_transpose()
		 */
		// VirtualRowvec<T_Scalar> transpose() const;

		/* TODO: use virtuals
		 * @copydoc standard_vector_docs::virtual_ctranspose()
		 */
		//VirtualRowvec<T_Scalar> ctranspose() const;

		/**
		 * @copydoc standard_vector_docs::virtual_conjugate()
		 */
		XxVector<T_Scalar> conjugate() const; // TODO: use virtuals

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
		::cla3p::Guard<XxVector<T_Scalar>> rblock(int_t ibgn, int_t ni) const;

		/**
		 * @copydoc standard_vector_docs::set_block()
		 */
		void setBlock(int_t ibgn, const XxVector<T_Scalar>& src);

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_XXVECTOR_HPP_