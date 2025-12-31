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

#ifndef CLA3P_TUPLE_HPP_
#define CLA3P_TUPLE_HPP_

/**
 * @file
 */

#include <ostream>

#include "cla3p/generic/coord.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

template <typename T_Int, typename T_Scalar>
class Tuple : public Coord<T_Int> {

	public:
		Tuple() { defaults(); }
		Tuple(T_Int r, T_Int c, T_Scalar v) : Coord<T_Int>(r, c) { setVal(v); }
		~Tuple() { clear(); }

		T_Scalar val() const { return m_val; }

		void clear()
		{
			Coord<T_Int>::clear();
			defaults(); 
		}

	private:
		T_Scalar m_val;

		void setVal(T_Scalar v) { m_val = v; }

		void defaults() { setVal(0); }
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::Tuple<T_Int,T_Scalar>& tuple)
{
  os << "(" << tuple.row() << ", " << tuple.col() << ", " << tuple.val() << ")";
  return os;
}

#endif // CLA3P_TUPLE_HPP_
