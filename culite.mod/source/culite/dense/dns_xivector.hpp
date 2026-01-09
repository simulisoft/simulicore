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
#include <cla3p/dense/dns_xxcontainer.hpp>
#include <cla3p/generic/guard.hpp>

#include "culite/types/integer.hpp"



/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The general purpose device dense vector class.
 */
template <typename T_Scalar>
class XiVector : public ::cla3p::Meta1D<int_t>, public ::cla3p::dns::XxContainer<T_Scalar> {

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
		 * @copydoc standard_docs::clear()
		 */
		void clear();

		/**
		 * @copydoc standard_docs::copy()
		 */
		XiVector<T_Scalar> copy() const;

		/**
		 * @copydoc standard_docs::rcopy()
		 */
		XiVector<T_Scalar> rcopy();

		/**
		 * @copydoc standard_docs::rcopy_const()
		 */
		::cla3p::Guard<XiVector<T_Scalar>> rcopy() const;

		/**
		 * @copydoc standard_docs::move()
		 */
		XiVector<T_Scalar> move();

		/**
		 * @copydoc standard_vector_docs::info()
		 */
		std::string info(const std::string& header = "") const;

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @copydoc standard_vector_docs::view()
		 */
		static ::cla3p::Guard<XiVector<T_Scalar>> view(int_t n, const T_Scalar *vals);

		/** @} */

	private:
		void moveFrom(XiVector<T_Scalar>& other);
		void copyFromExisting(const XiVector<T_Scalar>& other);
		void checker() const;

		virtual T_Scalar* allocateMemory(std::size_t numElements) override;
		virtual void freeMemory(T_Scalar *ptr) override;
};

/*-------------------------------------------------*/

// TODO: implement a XiVectorHost class that inherits from ::cla3p::dns::XiVector
//       and overrides allocators with pinned memory allocators

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

/* TODO: enable when ready
 * @ingroup culite_module_index_stream_operators
 * @brief Writes to os the contents of vec.
 */
//template <typename T_Scalar>
//std::ostream& operator<<(std::ostream& os, const culite::dns::XiVector<T_Scalar>& vec)
//{
//	vec.toStream(os);
//	return os;
//}

#endif // CULITE_DNS_XIVECTOR_HPP_
