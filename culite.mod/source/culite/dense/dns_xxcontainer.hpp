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

#ifndef CCULITE_DNS_XXCONTAINER_HPP_
#define CCULITE_DNS_XXCONTAINER_HPP_

/**
 * @file
 */

#include <iostream>
#include <cstddef>

#include <cla3p/dense/dns_xxcontainer_base.hpp>
#include <cla3p/generic/ownership.hpp>

#include "culite/support/imalloc.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The device dense container class.
 */
template <typename T_Scalar>
class XxContainer : public ::cla3p::dns::XxContainerBase<T_Scalar>, public ::cla3p::Ownership {

	public:
		using value_type = T_Scalar;

	protected:
		XxContainer() {}

		explicit XxContainer(std::size_t numElements)
			: ::cla3p::dns::XxContainerBase<T_Scalar>(numElements ? device_alloc_t<T_Scalar>(numElements) : nullptr),
			  Ownership(numElements ? true : false) {}

		explicit XxContainer(T_Scalar *vals, bool bind)
			: ::cla3p::dns::XxContainerBase<T_Scalar>(vals),
			  Ownership(vals ? bind : false) {}

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

	protected:
		void clear()
		{
			if(owner()) {
				device_free(this->values());
			} // owner
			::cla3p::dns::XxContainerBase<T_Scalar>::clear();
			Ownership::clear();
		}

	private:
		void moveFrom(XxContainer<T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				::cla3p::dns::XxContainerBase<T_Scalar>::operator=(std::move(other));
				Ownership::operator=(std::move(other));
				other.unbind();
				other.clear();
			} // do not apply on self
		}
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CULITE_DNS_XXCONTAINER_HPP_
