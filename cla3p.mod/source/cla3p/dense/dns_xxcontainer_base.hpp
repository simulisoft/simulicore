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

#ifndef CLA3P_DNS_XXCONTAINER_BASE_HPP_
#define CLA3P_DNS_XXCONTAINER_BASE_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense container base class.
 */
template <typename T_Scalar>
class XxContainerBase {

	public:
		using value_type = T_Scalar;

	protected:
		XxContainerBase()
		{
			defaults();
		}

		XxContainerBase(T_Scalar *vals)
			: XxContainerBase<T_Scalar>()
		{
			setValues(vals);
		}

		XxContainerBase(XxContainerBase<T_Scalar>&) = delete;
		XxContainerBase<T_Scalar>& operator=(XxContainerBase<T_Scalar>&) = delete;

		XxContainerBase(XxContainerBase<T_Scalar>&& other)
		{
			moveFrom(other);
		}

		XxContainerBase<T_Scalar>& operator=(XxContainerBase<T_Scalar>&& other)
		{
			moveFrom(other);
			return *this;
		}

		~XxContainerBase()
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
			defaults();
		}

	private:
		T_Scalar *m_values;

		void setValues(T_Scalar *vals) { m_values = vals; }

		void defaults() { setValues(nullptr); }

		void moveFrom(XxContainerBase<T_Scalar>& other)
		{
			if(this != &other) {
				clear();
				setValues(other.values());
				other.clear();
			} // do not apply on self
		}
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DNS_XXCONTAINER_BASE_HPP_
