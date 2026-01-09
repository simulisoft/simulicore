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

#ifndef CULITE_DNS_CXVECTOR_HPP_
#define CULITE_DNS_CXVECTOR_HPP_

/**
 * @file
 */

#include "culite/dense/dns_xxvector.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The device dense complex vector class.
 */
template <typename T_Scalar>
class CxVector : public XxVector<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		explicit CxVector(const XiVector<T_Scalar>& other);
		CxVector(XiVector<T_Scalar>&& other);

    	CxVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);
    	CxVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_vector_docs::constructor()
		 */
		CxVector();

		/**
		 * @copydoc standard_vector_docs::dim_constructor()
		 */
		explicit CxVector(int_t n);

		/**
		 * @copydoc standard_vector_docs::aux_constructor()
		 */
		explicit CxVector(int_t n, T_Scalar *vals, bool bind);

		/**
		 * @copydoc standard_docs::copy_constructor()
		 */
		CxVector(const CxVector<T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_constructor()
		 */
		CxVector(CxVector<T_Scalar>&& other) = default;

		/**
		 * @copydoc standard_vector_docs::destructor()
		 */
		~CxVector();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @copydoc standard_docs::copy_assignment()
		 */
		CxVector<T_Scalar>& operator=(const CxVector<T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_assignment()
		 */
		CxVector<T_Scalar>& operator=(CxVector<T_Scalar>&& other) = default;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/* // TODO: implement
		 * @brief Gets a copy of the real part of the vector.
		 * @return A copy of the real part of the vector.
		 */
		//XxVector<T_RScalar> real() const;

		/* // TODO: implement
		 * @brief Gets a copy of the imaginary part of the vector.
		 * @return A copy of the imaginary part of the vector.
		 */
		//XxVector<T_RScalar> imag() const;

		/** @} */

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_CXVECTOR_HPP_
