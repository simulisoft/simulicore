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
 * @details Represents a complex-valued dense vector stored on the GPU device.
 * @tparam T_Scalar The complex scalar type (e.g., std::complex<float>, std::complex<double>).
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
        CxVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

		CxVector(XiVector<T_Scalar>&& other);
    	CxVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);

        /**
         * @name Virtual Convertors
         * @{
         */

        template <typename T_Virtual>
        CxVector(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v) : XxVector<T_Scalar>(v) {}
        template <typename T_Virtual>
        CxVector<T_Scalar>& operator=(const alias::VirtualExpr_vec<T_Scalar,T_Virtual>& v) 
        { 
            XxVector<T_Scalar>::operator=(v); 
            return *this; 
        }

        /** @} */

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Creates an empty complex device vector with no allocated device memory.
		 */
		CxVector();

		/**
		 * @brief Dimension constructor.
		 * @details Creates a complex device vector of the specified size and allocates device memory.
		 * @param[in] n The number of elements in the vector.
		 */
		explicit CxVector(int_t n);

		/**
		 * @brief Auxiliary constructor.
		 * @details Creates a complex device vector using existing device memory.
		 * @param[in] n The number of elements in the vector.
		 * @param[in] vals Pointer to existing device memory.
		 * @param[in] bind If true, the vector does not take ownership of the memory.
		 */
		explicit CxVector(int_t n, T_Scalar *vals, bool bind);

		/**
		 * @brief Copy constructor.
		 * @details Creates a new complex device vector by copying another complex device vector.
		 * @param[in] other The complex device vector to copy.
		 */
		CxVector(const CxVector<T_Scalar>& other) = default;

		/**
		 * @brief Move constructor.
		 * @details Creates a new complex device vector by moving resources from another complex device vector.
		 * @param[in] other The complex device vector to move from.
		 */
		CxVector(CxVector<T_Scalar>&& other) = default;

		/**
		 * @brief Destructor.
		 * @details Destroys the complex device vector and releases allocated device memory.
		 */
		~CxVector();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another complex device vector to this complex device vector.
		 * @param[in] other The complex device vector to copy.
		 * @return Reference to this complex device vector.
		 */
		CxVector<T_Scalar>& operator=(const CxVector<T_Scalar>& other) = default;

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another complex device vector to this complex device vector.
		 * @param[in] other The complex device vector to move from.
		 * @return Reference to this complex device vector.
		 */
		CxVector<T_Scalar>& operator=(CxVector<T_Scalar>&& other) = default;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Extract the real part.
		 * @details Creates a new device vector containing the real part of each complex element.
		 * @return A device vector containing the real parts.
		 */
		XxVector<T_RScalar> real() const;

		/**
		 * @brief Extract the imaginary part.
		 * @details Creates a new device vector containing the imaginary part of each complex element.
		 * @return A device vector containing the imaginary parts.
		 */
		XxVector<T_RScalar> imag() const;	

		/** @} */

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_CXVECTOR_HPP_
