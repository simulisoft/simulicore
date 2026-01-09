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

#ifndef CLA3P_DNS_XXCONTAINER_HPP_
#define CLA3P_DNS_XXCONTAINER_HPP_

/**
 * @file
 */

#include <cstddef>

#include "cla3p/generic/ownership.hpp"
#include "cla3p/generic/guard.hpp"
#include "cla3p/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense container class.
 */
template <typename T_Scalar>
class XxContainer : public Ownership {

	public:
		using value_type = T_Scalar;

	protected:
		XxContainer()
		{
			defaults();
		}

		explicit XxContainer(std::size_t numElements)
			: XxContainer<T_Scalar>()
		{
			if(numElements) {
				T_Scalar *vals = allocateMemory(numElements);
				Ownership::operator=(Ownership(true));
				setValues(vals);
			}
		}

		explicit XxContainer(T_Scalar *vals, bool bind)
			: XxContainer<T_Scalar>()
		{
			if(vals) {
				Ownership::operator=(Ownership(bind));
				setValues(vals);
			}
		}

		XxContainer(XxContainer<T_Scalar>&) = delete;
		XxContainer<T_Scalar>& operator=(XxContainer<T_Scalar>&) = delete;

		XxContainer(XxContainer<T_Scalar>&& other)
		{
			moveFrom(other);
		}

		XxContainer<T_Scalar>& operator=(XxContainer<T_Scalar>&& other)
		{
			moveFrom(other);
			return *this;
		}

		~XxContainer()
		{
			clear();
		}

	public:

		/**
		 * @copydoc standard_docs::values()
		 */
		T_Scalar* values() { return m_values; }

		/**
		 * @copydoc standard_docs::values()
		 */
		const T_Scalar* values() const { return m_values; }

	protected:
		void clear()
		{
			if(owner()) {
				freeMemory(values());
			} // owner
			Ownership::clear();
			defaults();
		}

	private:
		T_Scalar *m_values;

		void setValues(T_Scalar *vals) { m_values = vals; }

		void defaults() { setValues(nullptr); }

		void moveFrom(XxContainer<T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				Ownership::operator=(std::move(other));
				setValues(other.values());
				other.unbind();
				other.clear();
			} // do not apply on self
		}

		virtual T_Scalar* allocateMemory(std::size_t numElements)
		{
			return i_malloc<T_Scalar>(numElements);
		}

		virtual void freeMemory(T_Scalar *ptr)
		{
			i_free(ptr);
		}
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DNS_XXCONTAINER_HPP_
