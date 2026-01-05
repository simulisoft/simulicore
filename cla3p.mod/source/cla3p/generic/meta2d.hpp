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

#ifndef CLA3P_META2D_HPP_
#define CLA3P_META2D_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Two-dimensional metadata container.
 * @tparam T_Int Integer type for storing dimension values.
 * @details Provides storage and accessors for matrix dimensions (rows and columns).
 *          This class serves as a base for dimension metadata in dense matrix structures.
 */
template <typename T_Int>
class Meta2D {

	public:
		Meta2D() { defaults(); }
		Meta2D(T_Int m, T_Int n)
		{
			setNrows(m);
			setNcols(n);
		}
		~Meta2D() { clear(); }

		/**
		 * @brief Number of rows.
		 * @details Returns the number of rows stored in @p (*this).
		 * @return The row dimension.
		 */
		virtual T_Int nrows() const { return m_nrows; }

		/**
		 * @brief Number of columns.
		 * @details Returns the number of columns stored in @p (*this).
		 * @return The column dimension.
		 */
		virtual T_Int ncols() const { return m_ncols; }

		/**
		 * @brief Tests whether the object has zero dimensions.
		 * @details Returns @c true if either dimension is zero.
		 * @return @c true if the object is empty, @c false otherwise.
		 */
		bool empty() const { return !(nrows() && ncols()); }

		/**
		 * @brief Boolean conversion operator.
		 * @details Enables implicit conversion to @c bool; returns @c true if object
		 *          has non-zero dimensions.
		 * @return @c true if the object is not empty, @c false otherwise.
		 */
		operator bool() const { return !empty(); }

	protected:
		void clear() { defaults(); }

	private:
		T_Int m_nrows;
		T_Int m_ncols;

		void setNrows(T_Int m) { m_nrows = m; }
		void setNcols(T_Int n) { m_ncols = n; }

		void defaults()
		{
			setNrows(0);
			setNcols(0);
		}
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_META2D_HPP_
