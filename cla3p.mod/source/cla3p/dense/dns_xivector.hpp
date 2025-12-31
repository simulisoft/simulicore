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

#ifndef CLA3P_DNS_XIVECTOR_HPP_
#define CLA3P_DNS_XIVECTOR_HPP_

/**
 * @file
 */

#include <string>
#include <ostream>

#include "cla3p/generic/guard.hpp"
#include "cla3p/generic/meta1d.hpp"
#include "cla3p/dense/dns_xxcontainer.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

namespace prm { template <typename T_Int> class PxMatrix; }

/*-------------------------------------------------*/
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The general purpose dense vector class.
 */
template <typename T_Scalar>
class XiVector : public Meta1D<int_t>, public XxContainer<T_Scalar> {

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
		 * @name Operators
		 * @{
		 */

		/**
		 * @copydoc standard_docs::index_operator_1d()
		 */
		T_Scalar& operator()(int_t i);

		/**
		 * @copydoc standard_docs::index_operator_1d()
		 */
		const T_Scalar& operator()(int_t i) const;

		/** @} */

		/** 
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @copydoc standard_docs::clear()
		 */
		void clear();

		/**
		 * @copydoc standard_vector_docs::fill()
		 */
		void fill(T_Scalar val);

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
		Guard<XiVector<T_Scalar>> rcopy() const;

		/**
		 * @copydoc standard_docs::move()
		 */
		XiVector<T_Scalar> move();

		/**
		 * @copydoc standard_vector_docs::info()
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @copydoc standard_docs::print()
		 */
		void print(std::streamsize prec = 0) const;

		/**
		 * @copydoc standard_docs::toStream()
		 */
		void toStream(std::ostream& os, std::streamsize prec = 0) const;

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @copydoc standard_vector_docs::view()
		 */
		static Guard<XiVector<T_Scalar>> view(int_t n, const T_Scalar *vals);

		/** @} */

	private:
		void moveFrom(XiVector<T_Scalar>& other);
		void copyFromExisting(const XiVector<T_Scalar>& other);
		void checker() const;
};

/*-------------------------------------------------*/
} // namespace dns
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of vec.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::dns::XiVector<T_Scalar>& vec)
{
	vec.toStream(os);
	return os;
}

#endif // CLA3P_DNS_XIVECTOR_HPP_
