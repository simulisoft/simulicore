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
        XxVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

		XxVector(XiVector<T_Scalar>&& other);
		XxVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Creates an empty device vector with no allocated device memory.
		 */
		XxVector();

		/**
		 * @brief Dimension constructor.
		 * @details Creates a device vector of the specified size and allocates device memory.
		 * @param[in] n The number of elements in the vector.
		 */
		explicit XxVector(int_t n);

		/**
		 * @brief Auxiliary constructor.
		 * @details Creates a device vector using existing device memory.
		 * @param[in] n The number of elements in the vector.
		 * @param[in] vals Pointer to existing device memory.
		 * @param[in] bind If true, the vector does not take ownership of the memory.
		 */
		explicit XxVector(int_t n, T_Scalar *vals, bool bind);

		/**
		 * @brief Copy constructor.
		 * @details Creates a new device vector by copying another device vector.
		 * @param[in] other The device vector to copy.
		 */
		XxVector(const XxVector<T_Scalar>& other) = default;

		/**
		 * @brief Move constructor.
		 * @details Creates a new device vector by moving resources from another device vector.
		 * @param[in] other The device vector to move from.
		 */
		XxVector(XxVector<T_Scalar>&& other) = default;

		/**
		 * @brief Destructor.
		 * @details Destroys the device vector and releases allocated device memory.
		 */
		~XxVector();

		/** @} */

		/** 
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another device vector to this device vector.
		 * @param[in] other The device vector to copy.
		 * @return Reference to this device vector.
		 */
		XxVector<T_Scalar>& operator=(const XxVector<T_Scalar>& other) = default;

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another device vector to this device vector.
		 * @param[in] other The device vector to move from.
		 * @return Reference to this device vector.
		 */
		XxVector<T_Scalar>& operator=(XxVector<T_Scalar>&& other) = default;

		/**
		 * @brief Unary negation operator.
		 * @details Returns a negated copy of the device vector.
		 * @return A device vector containing the negated elements.
		 */
		XxVector<T_Scalar> operator-() const; // TODO: use virtuals

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Scale the device vector in-place.
		 * @details Multiplies all elements of the device vector by a scalar value.
		 * @param[in] val The scalar value to multiply by.
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
		 * @brief Compute the complex conjugate.
		 * @details Returns a device vector containing the complex conjugate of each element.
		 * @return A device vector with conjugated elements.
		 */
		XxVector<T_Scalar> conjugate() const; // TODO: use virtuals

		/**
		 * @brief Conjugate the device vector in-place.
		 * @details Replaces all elements with their complex conjugates.
		 */
		void iconjugate();

		/**
		 * @brief Compute the 1-norm.
		 * @details Computes the sum of absolute values of all elements.
		 * @return The 1-norm of the device vector.
		 */
		T_RScalar normOne() const;

		/**
		 * @brief Compute the infinity norm.
		 * @details Computes the maximum absolute value of all elements.
		 * @return The infinity norm of the device vector.
		 */
		T_RScalar normInf() const;

		/**
		 * @brief Compute the Euclidean norm (2-norm).
		 * @details Computes the square root of the sum of squared absolute values.
		 * @return The Euclidean norm of the device vector.
		 */
		T_RScalar normEuc() const;

		/**
		 * @brief Extract a block as a new device vector.
		 * @details Creates a new device vector containing a deep copy of a contiguous block of elements.
		 * @param[in] ibgn The starting index of the block.
		 * @param[in] ni The number of elements in the block.
		 * @return A new device vector containing the block's data.
		 */
		XxVector<T_Scalar> block(int_t ibgn, int_t ni) const;

		/**
		 * @brief Extract a reference block.
		 * @details Creates a device vector that references a contiguous block of this vector's device memory.
		 * @param[in] ibgn The starting index of the block.
		 * @param[in] ni The number of elements in the block.
		 * @return A device vector that shares device memory with this vector.
		 */
		XxVector<T_Scalar> rblock(int_t ibgn, int_t ni);

		/**
		 * @brief Extract a guarded reference block.
		 * @details Creates a guarded device vector that references a contiguous block of this vector's device memory.
		 * @param[in] ibgn The starting index of the block.
		 * @param[in] ni The number of elements in the block.
		 * @return A guarded device vector that shares device memory with this vector.
		 */
		::cla3p::Guard<XxVector<T_Scalar>> rblock(int_t ibgn, int_t ni) const;

		/**
		 * @brief Set a block of elements.
		 * @details Copies elements from the source device vector to a contiguous block of this vector.
		 * @param[in] ibgn The starting index where the block will be written.
		 * @param[in] src The source device vector to copy from.
		 */
		void setBlock(int_t ibgn, const XxVector<T_Scalar>& src);

};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_XXVECTOR_HPP_