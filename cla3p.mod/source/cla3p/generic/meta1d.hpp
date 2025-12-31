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

#ifndef CLA3P_META1D_HPP_
#define CLA3P_META1D_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The single dimension metadata class.
 */
template <typename T_Int>
class Meta1D {

	public:
		Meta1D() { defaults(); }
		Meta1D(T_Int n) { setSize(n); }
		~Meta1D() { clear(); }

		/**
		 * @brief The vector size.
		 * @return The number of entries in `(*this)`.
		 */
		T_Int size() const { return m_size; }

		/**
		 * @brief Test whether object is empty.
		 * @return true if object dimension is zero, false otherwise
		 */
		bool empty() const { return !size(); }

		/**
		 * @copydoc empty()
		 */
		bool operator!() const { return empty(); } 

		/**
		 * @brief Test whether object is empty.
		 * @return true if object is not empty, false otherwise
		 */
		operator bool() const { return !empty(); }

	protected:
		void clear() { defaults(); }

	private:
		T_Int m_size;

		void defaults() { setSize(0); }

		void setSize(T_Int n) { m_size = n; }
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_META1D_HPP_
