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

#ifndef CLA3P_COORD_HPP_
#define CLA3P_COORD_HPP_

/**
 * @file
 */

#include <ostream>

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Two-dimensional coordinate container.
 * @tparam T_Int Integer type for coordinate indices.
 * @details Stores a row and column index pair representing a position
 *          within a two-dimensional structure.
 */
template <typename T_Int>
class Coord {

    public:
        /**
         * @brief Default constructor.
         * @details Constructs a coordinate initialized to (0, 0).
         */
        Coord() { defaults(); }

        /**
         * @brief Coordinate constructor.
         * @details Constructs a coordinate with specified row and column indices.
         * @param[in] r The row index.
         * @param[in] c The column index.
         */
        Coord(T_Int r, T_Int c)
        {
            setRow(r);
            setCol(c);
        }

        /**
         * @brief Destructor.
         * @details Destroys the coordinate object.
         */
        ~Coord() { clear(); }

        /**
         * @brief Row index accessor.
         * @details Returns the row coordinate value.
         * @return The row index.
         */
        T_Int row() const { return m_row; }

        /**
         * @brief Column index accessor.
         * @details Returns the column coordinate value.
         * @return The column index.
         */
        T_Int col() const { return m_col; }

    protected:
        void clear() { defaults(); }

    private:
        T_Int m_row;
        T_Int m_col;

        void setRow(T_Int r) { m_row = r; }
        void setCol(T_Int c) { m_col = c; }

        void defaults()
        {
            setRow(0);
            setCol(0);
        }
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @brief Stream insertion operator for Coord.
 * @tparam T_Int Integer type for coordinate indices.
 * @details Outputs the coordinate in the format (row, col) to the stream.
 * @param[in,out] os The output stream.
 * @param[in] coord The coordinate to output.
 * @return Reference to the output stream.
 */
template <typename T_Int>
std::ostream& operator<<(std::ostream& os, const cla3p::Coord<T_Int>& coord)
{
  os << "(" << coord.row() << ", " << coord.col() << ")";
  return os;
}

#endif // CLA3P_COORD_HPP_
