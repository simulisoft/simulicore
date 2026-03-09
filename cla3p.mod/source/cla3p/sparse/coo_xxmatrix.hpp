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

#ifndef CLA3P_COO_XXMATRIX_HPP_
#define CLA3P_COO_XXMATRIX_HPP_

/**
 * @file
 */

#include <ostream>
#include <string>
#include <vector>

#include "cla3p/types/enums.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/generic/tuple.hpp"
#include "cla3p/sparse/csr_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace coo {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse matrix class (coordinate format).
 */
template <typename T_Int, typename T_Scalar>
class XxMatrix : public MatrixMeta<T_Int> {

    private:
        using TupleVec = std::vector<Tuple<T_Int,T_Scalar>>;

    public:
        using index_type = T_Int;
        using value_type = T_Scalar;

    public:

        /**
         * @name Constructors
         * @{
         */

        /**
         * @brief Default constructor.
         * @details Creates an empty sparse matrix with no allocated memory.
         */
        XxMatrix();

        /**
         * @brief Dimension constructor.
         * @details Creates a sparse matrix of the specified dimensions in coordinate format.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] pr The matrix property (default: General).
         */
        explicit XxMatrix(T_Int nr, T_Int nc, const Property& pr = Property::General());

        /**
         * @brief Copy constructor.
         * @details Creates a new sparse matrix by copying another sparse matrix.
         * @param[in] other The sparse matrix to copy.
         */
        XxMatrix(const XxMatrix<T_Int,T_Scalar>& other) = default;

        /**
         * @brief Move constructor.
         * @details Creates a new sparse matrix by moving resources from another sparse matrix.
         * @param[in] other The sparse matrix to move from.
         */
        XxMatrix(XxMatrix<T_Int,T_Scalar>&& other) = default;

        /**
         * @brief Destructor.
         * @details Destroys the sparse matrix and releases allocated memory.
         */
        ~XxMatrix();

        /** @} */

        /** 
         * @name Operators
         * @{
         */

        /**
         * @brief Copy assignment operator.
         * @details Copies the contents of another sparse matrix to this sparse matrix.
         * @param[in] other The sparse matrix to copy.
         * @return Reference to this sparse matrix.
         */
        XxMatrix<T_Int,T_Scalar>& operator=(const XxMatrix<T_Int,T_Scalar>& other) = default;

        /**
         * @brief Move assignment operator.
         * @details Moves resources from another sparse matrix to this sparse matrix.
         * @param[in] other The sparse matrix to move from.
         * @return Reference to this sparse matrix.
         */
        XxMatrix<T_Int,T_Scalar>& operator=(XxMatrix<T_Int,T_Scalar>&& other) = default;

        /** @} */

        /** 
         * @name Arguments
         * @{
         */

        /**
         * @brief Get the number of non-zero elements.
         * @details Returns the number of stored non-zero elements in the sparse matrix.
         * @return The number of non-zero elements.
         */
        T_Int nnz() const;

        /** @} */

        /** 
         * @name Public Member Functions
         * @{
         */

        /**
         * @brief Clear the sparse matrix.
         * @details Releases all memory and resets the matrix to an empty state.
         */
        void clear();

        /**
         * @brief Reserve storage for non-zero elements.
         * @details Preallocates memory for the specified number of non-zero elements to avoid reallocations.
         * @param[in] nz The number of non-zero elements to reserve space for.
         */
        void reserve(T_Int nz);

        /**
         * @brief Insert a non-zero element using a tuple.
         * @details Inserts a non-zero element specified as a tuple (row, column, value) into the sparse matrix.
         * @param[in] tuple The tuple containing the row index, column index, and value.
         */
        void insert(const Tuple<T_Int,T_Scalar>& tuple);

        /**
         * @brief Insert a non-zero element.
         * @details Inserts a non-zero element at the specified row and column with the given value.
         * @param[in] i The row index.
         * @param[in] j The column index.
         * @param[in] v The value to insert.
         */
        void insert(T_Int i, T_Int j, T_Scalar v);

        /**
         * @brief Get information about the sparse matrix.
         * @details Returns a string containing information about the matrix's dimensions, non-zeros, and properties.
         * @param[in] header Optional header string to prepend to the information.
         * @return A string containing matrix information.
         */
        std::string info(const std::string& header = "") const;

        /**
         * @brief Output the sparse matrix to a stream.
         * @details Writes the sparse matrix elements to the specified output stream with specified precision.
         * @param[in,out] os The output stream to write to.
         * @param[in] prec The output precision (default: 0 for default precision).
         */
        void toStream(std::ostream& os, std::streamsize prec = 0) const;

        /**
         * @brief Convert to compressed sparse column (CSC) format.
         * @details Converts the coordinate format sparse matrix to CSC format.
         * @param[in] duplicatePolicy The policy for handling duplicate entries (default: Sum).
         * @return A sparse matrix in CSC format.
         */
        csc::XxMatrix<T_Int,T_Scalar> toCsc(dup_t duplicatePolicy = dup_t::Sum) const;

        /**
         * @brief Convert to compressed sparse row (CSR) format.
         * @details Converts the coordinate format sparse matrix to CSR format.
         * @param[in] duplicatePolicy The policy for handling duplicate entries (default: Sum).
         * @return A sparse matrix in CSR format.
         */
        csr::XxMatrix<T_Int,T_Scalar> toCsr(dup_t duplicatePolicy = dup_t::Sum) const;

        /** @} */

    private:
        TupleVec m_tuples;

        TupleVec& tupleVec();
        const TupleVec& tupleVec() const;

        void checker() const;
};

/*-------------------------------------------------*/
} // namespace coo
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 */
template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::coo::XxMatrix<T_Int,T_Scalar>& mat)
{
    mat.toStream(os);
    return os;
}

#endif // CLA3P_COO_XXMATRIX_HPP_
