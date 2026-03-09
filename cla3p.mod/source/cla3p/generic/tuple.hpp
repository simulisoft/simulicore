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

/**
 * @nosubgrouping
 * @brief Coordinate-value triple container.
 * @tparam T_Int Integer type for coordinate indices.
 * @tparam T_Scalar Scalar type for the stored value.
 * @details Extends @ref Coord with an associated scalar value, representing
 *          a position and its corresponding data within a two-dimensional structure.
 */
template <typename T_Int, typename T_Scalar>
class Tuple : public Coord<T_Int> {

    public:
        /**
         * @brief Default constructor.
         * @details Constructs a tuple initialized to (0, 0, 0).
         */
        Tuple() { defaults(); }

        /**
         * @brief Value constructor.
         * @details Constructs a tuple with specified row, column, and value.
         * @param[in] r The row index.
         * @param[in] c The column index.
         * @param[in] v The scalar value.
         */
        Tuple(T_Int r, T_Int c, T_Scalar v) : Coord<T_Int>(r, c) { setVal(v); }

        /**
         * @brief Destructor.
         * @details Destroys the tuple object.
         */
        ~Tuple() { clear(); }

        /**
         * @brief Value accessor.
         * @details Returns the scalar value stored in the tuple.
         * @return The tuple value.
         */
        T_Scalar val() const { return m_val; }

        /**
         * @brief Clears the tuple.
         * @details Resets the coordinate and value to their default states (0, 0, 0).
         */
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

/**
 * @brief Stream insertion operator for Tuple.
 * @tparam T_Int Integer type for coordinate indices.
 * @tparam T_Scalar Scalar type for the stored value.
 * @details Outputs the tuple in the format (row, col, val) to the stream.
 * @param[in,out] os The output stream.
 * @param[in] tuple The tuple to output.
 * @return Reference to the output stream.
 */
template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::Tuple<T_Int,T_Scalar>& tuple)
{
  os << "(" << tuple.row() << ", " << tuple.col() << ", " << tuple.val() << ")";
  return os;
}

#endif // CLA3P_TUPLE_HPP_
