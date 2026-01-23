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
 * @details Represents a complex-valued dense matrix stored on the GPU device in column-major format.
 * @tparam T_Scalar The complex scalar type (e.g., std::complex<float>, std::complex<double>).
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
		 * @brief Default constructor.
		 * @details Creates an empty device complex matrix with no allocated device memory.
		 */
		CxMatrix();

		/**
		 * @brief Dimension constructor.
		 * @details Creates a device complex matrix of the specified dimensions and allocates device memory.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] pr The matrix property (default: General).
		 */
		explicit CxMatrix(int_t nr, int_t nc, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @brief Auxiliary constructor.
		 * @details Creates a device complex matrix using existing device memory.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] vals Pointer to existing device memory.
		 * @param[in] ldv The leading dimension.
		 * @param[in] bind If true, the matrix does not take ownership of the memory.
		 * @param[in] pr The matrix property (default: General).
		 */
		explicit CxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @brief Copy constructor.
		 * @details Creates a new device complex matrix by copying another device complex matrix.
		 * @param[in] other The device complex matrix to copy.
		 */
		CxMatrix(const CxMatrix<T_Scalar>& other) = default;

		/**
		 * @brief Move constructor.
		 * @details Creates a new device complex matrix by moving resources from another device complex matrix.
		 * @param[in] other The device complex matrix to move from.
		 */
		CxMatrix(CxMatrix<T_Scalar>&& other) = default;

		/**
		 * @brief Destructor.
		 * @details Destroys the device complex matrix and releases allocated device memory.
		 */
		~CxMatrix();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another device complex matrix to this device complex matrix.
		 * @param[in] other The device complex matrix to copy.
		 * @return Reference to this device complex matrix.
		 */
		CxMatrix<T_Scalar>& operator=(const CxMatrix<T_Scalar>& other) = default;

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another device complex matrix to this device complex matrix.
		 * @param[in] other The device complex matrix to move from.
		 * @return Reference to this device complex matrix.
		 */
		CxMatrix<T_Scalar>& operator=(CxMatrix<T_Scalar>&& other) = default;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Extract the real part.
		 * @details Creates a new device matrix containing the real parts of all complex elements.
		 * @return A device matrix containing the real components.
		 */
		XxMatrix<T_RScalar> real() const;

		/**
		 * @brief Extract the imaginary part.
		 * @details Creates a new device matrix containing the imaginary parts of all complex elements.
		 * @return A device matrix containing the imaginary components.
		 */
		XxMatrix<T_RScalar> imag() const;

		/** @} */

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_CXMATRIX_HPP_