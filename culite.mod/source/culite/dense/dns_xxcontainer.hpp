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

#ifndef CULITE_DNS_XXCONTAINER_HPP_
#define CULITE_DNS_XXCONTAINER_HPP_

/**
 * @file
 */

#include <cstddef>
#include <cla3p/generic/ownership.hpp>

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense container class.
 */
template <typename T_Scalar>
class XxContainer : public ::cla3p::Ownership {

	public:
		using value_type = T_Scalar;

	protected:
		XxContainer();
		explicit XxContainer(std::size_t numElements);
		explicit XxContainer(T_Scalar *vals, bool bind);

		XxContainer(XxContainer<T_Scalar>&) = delete;
		XxContainer<T_Scalar>& operator=(XxContainer<T_Scalar>&) = delete;

		XxContainer(XxContainer<T_Scalar>&&);
		XxContainer<T_Scalar>& operator=(XxContainer<T_Scalar>&&);

		~XxContainer();

	public:

		/**
		 * @copydoc standard_docs::values()
		 */
		T_Scalar* values();

		/**
		 * @copydoc standard_docs::values()
		 */
		const T_Scalar* values() const;

	protected:
		void clear();

	private:
		T_Scalar *m_values;

		void setValues(T_Scalar *vals);

		void defaults();

		void moveFrom(XxContainer<T_Scalar>& other);
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DNS_XXCONTAINER_HPP_
