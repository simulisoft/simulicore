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

#ifndef CULITE_GUARD_HPP_
#define CULITE_GUARD_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/
/**
* @brief Encloses an object and prevents its contents from being changed.
* @details Some functions need to return objects that contain immutable data.
*          The contents of these objects must be protected.
*          The Guard class is used for preventing those data from being exposed to change.
*/

template <typename T_Object>
class Guard {

	public:

		/**
		 * @copybrief standard_docs::constructor()
		 * @details Constructs an empty guard.
		 */
		Guard() {}

		/**
		 * @brief Destroys the guard.
		 * @details Detaches contents and destroys guard.
		 */
		~Guard() { clear(); }

		/**
		 * @copybrief standard_docs::copy_constructor()
		 * @details Constructs a guard with a referenced copy of the contents of `other`.
		 */
		Guard(const Guard<T_Object>& other) { setLocalObj(other.get()); }

		/**
		 * @copybrief standard_docs::copy_assignment()
		 * @details Replaces the contents of guard with a referenced copy of the contents of `other`.
		 */
		Guard<T_Object>& operator=(const Guard<T_Object>& other)
		{
			setLocalObj(other.get());
			return *this;
		}

		/**
		 * @copybrief standard_docs::move_convertor()
		 */
		template <typename U_Object>
		Guard(Guard<U_Object>&& other) 
			: Guard(static_cast<const T_Object&>(other.get())) {} 

		/**
		 * @brief The input constructor.
		 * @details Constructs guard with a referenced copy of `obj`.
		 */
		explicit Guard(const T_Object& obj) { setLocalObj(obj); }

		/**
		 * @brief Clears the guard.
		 * @details Detaches contents and resets guard.
		 */
		void clear() { m_obj.clear(); }

		/**
		 * @brief The guarded object.
		 * @return A constant reference to the object being guarded.
		 */
		const T_Object& get() const { return m_obj; }

	private:
		T_Object m_obj;

		void setLocalObj(const T_Object& obj)
		{
			m_obj.clear();
			m_obj = const_cast<T_Object&>(obj).rcopy(); 
		}
};
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_GUARD_HPP_
