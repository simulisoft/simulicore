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

#ifndef CLA3P_MATRIX_META_HPP_
#define CLA3P_MATRIX_META_HPP_

/**
 * @file
 */

#include "cla3p/generic/meta2d.hpp"
#include "cla3p/types/property.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The matrix metadata class.
 */
template <typename T_Int>
class MatrixMeta : public Meta2D<T_Int> {

	public:
		MatrixMeta() {}
		MatrixMeta(T_Int nr, T_Int nc, const Property& pr) : Meta2D<T_Int>(nr, nc) { setProp(pr); }
		~MatrixMeta() { clear(); }

		/**
		 * @brief The matrix property.
		 * @details Gets the defined property for the matrix.
		 */
		const Property& prop() const { return m_prop; }

	protected:
		void clear()
		{
			Meta2D<T_Int>::clear();
			m_prop.clear();
		}

		void setProp(const Property& pr) { m_prop = pr; }

	private:
		Property m_prop;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_MATRIX_META_HPP_
