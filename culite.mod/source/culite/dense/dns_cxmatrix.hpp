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

#ifndef CULITE_DNS_CXMATRIX_HPP_
#define CULITE_DNS_CXMATRIX_HPP_

/**
 * @file
 */

#include "culite/dense/dns_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The device dense complex matrix class.
 */
template <typename T_Scalar>
class CxMatrix : public XxMatrix<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		explicit CxMatrix(const XxMatrix<T_Scalar>& other);
		CxMatrix(XxMatrix<T_Scalar>&& other);

		CxMatrix<T_Scalar>& operator=(const XxMatrix<T_Scalar>& other);
		CxMatrix<T_Scalar>& operator=(XxMatrix<T_Scalar>&& other);

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @copydoc standard_matrix_docs::constructor()
		 */
		CxMatrix();

		/**
		 * @copydoc standard_matrix_docs::dim_constructor()
		 */
		explicit CxMatrix(int_t nr, int_t nc, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @copydoc standard_matrix_docs::aux_constructor()
		 */
		explicit CxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @copydoc standard_docs::copy_constructor()
		 */
		CxMatrix(const CxMatrix<T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_constructor()
		 */
		CxMatrix(CxMatrix<T_Scalar>&& other) = default;

		/**
		 * @copydoc standard_matrix_docs::destructor()
		 */
		~CxMatrix();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @copydoc standard_docs::copy_assignment()
		 */
		CxMatrix<T_Scalar>& operator=(const CxMatrix<T_Scalar>& other) = default;

		/**
		 * @copydoc standard_docs::move_assignment()
		 */
		CxMatrix<T_Scalar>& operator=(CxMatrix<T_Scalar>&& other) = default;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Gets a copy of the real part of the matrix.
		 * @return A copy of the real part of the matrix.
		 */
		XxMatrix<T_RScalar> real() const;

		/**
		 * @brief Gets a copy of the imaginary part of the matrix.
		 * @return A copy of the imaginary part of the matrix.
		 */
		XxMatrix<T_RScalar> imag() const;

		/** @} */

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_CXMATRIX_HPP_