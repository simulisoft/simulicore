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
 * @brief One-dimensional metadata container.
 * @tparam T_Int Integer type for storing the dimension value.
 * @details Provides storage and accessors for vector size.
 *          This class serves as a base for dimension metadata in dense vector structures.
 */
template <typename T_Int>
class Meta1D {

	public:
		Meta1D() { defaults(); }
		Meta1D(T_Int n) { setSize(n); }
		~Meta1D() { clear(); }

		/**
		 * @brief Vector size.
		 * @details Returns the number of entries stored in @p (*this).
		 * @return The size dimension.
		 */
		T_Int size() const { return m_size; }

		/**
		 * @brief Tests whether the object has zero dimension.
		 * @details Returns @c true if the size is zero.
		 * @return @c true if the object is empty, @c false otherwise.
		 */
		bool empty() const { return !size(); }

		/* No need for impl
		 * @brief Logical negation operator.
		 * @details Equivalent to @ref empty(); returns @c true if object has zero dimension.
		 * @return @c true if the object is empty, @c false otherwise.
		 */
		// bool operator!() const { return empty(); } 

		/**
		 * @brief Boolean conversion operator.
		 * @details Enables implicit conversion to @c bool; returns @c true if object
		 *          has non-zero dimension.
		 * @return @c true if the object is not empty, @c false otherwise.
		 */
		explicit operator bool() const { return !empty(); }

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
