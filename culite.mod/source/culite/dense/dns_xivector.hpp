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

#ifndef CULITE_DNS_XIVECTOR_HPP_
#define CULITE_DNS_XIVECTOR_HPP_

/**
 * @file
 */

#include <string>
#include <ostream>
#include <cla3p/generic/meta1d.hpp>
#include <cla3p/generic/guard.hpp>
#include <cla3p/dense/dns_xivector.hpp>

#include "culite/types/integer.hpp"
#include "culite/dense/dns_xxcontainer.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The general purpose device dense vector class.
 * @details Represents a dense vector stored on the GPU device.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 */
template <typename T_Scalar>
class XiVector : public ::cla3p::Meta1D<int_t>, public XxContainer<T_Scalar> {

	private:
		using T_Cla3pScalar = typename TypeTraits<T_Scalar>::cla3p_type;

	public:
		XiVector();
		explicit XiVector(int_t n);
		explicit XiVector(int_t n, T_Scalar *vals, bool bind);
		~XiVector();

		XiVector(const XiVector<T_Scalar>& other);
		XiVector<T_Scalar>& operator=(const XiVector<T_Scalar>& other);

		XiVector(XiVector<T_Scalar>&& other);
		XiVector<T_Scalar>& operator=(XiVector<T_Scalar>&& other);

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Clear the device vector.
		 * @details Releases all device memory and resets the vector to an empty state.
		 */
		void clear();

        /**
         * @brief Fill all elements with a value.
         * @details Sets all elements in the vector to the specified scalar value.
         * @param[in] val The scalar value to fill with.
         */
        void fill(T_Scalar val);

		/**
		 * @brief Create a deep copy of the device vector.
		 * @details Creates a new device vector with its own device memory allocation
		 *          and copies all elements from this vector to the new vector.
		 * @return A new device vector containing a copy of this vector's data.
		 */
		XiVector<T_Scalar> copy() const;

		/**
		 * @brief Create a reference copy (shallow copy) of the device vector.
		 * @details Creates a new device vector object that references the same device memory
		 *          as this vector. Changes to either vector will affect both.
		 * @return A device vector that shares device memory with this vector.
		 */
		XiVector<T_Scalar> rcopy();

		/**
		 * @brief Create a guarded reference copy (shallow copy) of the device vector.
		 * @details Creates a guarded device vector object that references the same device memory
		 *          as this vector. The guard ensures the reference is read-only.
		 * @return A guarded device vector that shares device memory with this vector.
		 */
		::cla3p::Guard<XiVector<T_Scalar>> rcopy() const;

		/**
		 * @brief Move the device vector's resources.
		 * @details Transfers ownership of the device memory to a new vector object,
		 *          leaving this vector in an empty state.
		 * @return A new device vector containing this vector's device memory.
		 */
		XiVector<T_Scalar> move();

		/**
		 * @brief Get information about the device vector.
		 * @details Returns a string containing information about the vector's size and properties.
		 * @param[in] header Optional header string to prepend to the information.
		 * @return A string containing vector information.
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @brief Copies the device vector to a host vector.
		 * @details Performs a device-to-host copy operation.
		 *          If @p dest is empty, it is resized to match the size of the device vector.
		 *          Otherwise, @p dest must be compatible with the device vector.
		 * @param[out] dest The destination host vector.
		 */
		void copyToHost(::cla3p::dns::XiVector<T_Cla3pScalar>& dest) const;

		/**
		 * @brief Copies a host vector to the device vector.
		 * @details Performs a host-to-device copy operation.
		 *          If the device vector is empty, it is resized to match the size of @p src.
		 *          Otherwise, the device vector must be compatible with @p src.
		 * @param[in] src The source host vector.
		 */
		void copyFromHost(const ::cla3p::dns::XiVector<T_Cla3pScalar>& src);

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @brief Create a view of existing device memory.
		 * @details Creates a guarded device vector that references existing device memory
		 *          without taking ownership. The memory must remain valid for the lifetime
		 *          of the returned view.
		 * @param[in] n The number of elements in the vector.
		 * @param[in] vals Pointer to the device memory.
		 * @return A guarded device vector that views the specified device memory.
		 */
		static ::cla3p::Guard<XiVector<T_Scalar>> view(int_t n, const T_Scalar *vals);

		/** @} */

	private:
		XiVector<T_Scalar>& moveFrom(XiVector<T_Scalar>& other);
		XiVector<T_Scalar>& copyFromExisting(const XiVector<T_Scalar>& other);
		void checker() const;
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_device2host_operators
 * @brief Copies the contents of src to dest.
 * @details Performs a device-to-host copy from @p src to @p dest.
 *          If @p dest is empty, it is resized to match @p src, otherwise it must be compatible with @p src.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] src The source device vector.
 * @param[out] dest The destination host vector.
 */
template <typename T_Scalar>
void operator>>(const culite::dns::XiVector<T_Scalar>& src,
	            ::cla3p::dns::XiVector<typename culite::TypeTraits<T_Scalar>::cla3p_type>& dest)
{
	src.copyToHost(dest);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_host2device_operators
 * @brief Copies the contents of src to dest.
 * @details Performs a host-to-device copy from @p src to @p dest.
 *          If @p dest is empty, it is resized to match @p src, otherwise it must be compatible with @p src.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] src The source host vector.
 * @param[out] dest The destination device vector.
 */
template <typename T_Scalar>
void operator>>(const ::cla3p::dns::XiVector<typename culite::TypeTraits<T_Scalar>::cla3p_type>& src,
	            culite::dns::XiVector<T_Scalar>& dest)
{
	dest.copyFromHost(src);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Writes to os the contents of vec.
 * @details Copies the device vector to host and writes it to the output stream.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in,out] os The output stream.
 * @param[in] vec The device vector to output.
 * @return The output stream.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const culite::dns::XiVector<T_Scalar>& vec)
{
	::cla3p::dns::XiVector<typename culite::TypeTraits<T_Scalar>::cla3p_type> hostVec;
	vec >> hostVec;
	os << hostVec;
	return os;
}

/*-------------------------------------------------*/

#endif // CULITE_DNS_XIVECTOR_HPP_
