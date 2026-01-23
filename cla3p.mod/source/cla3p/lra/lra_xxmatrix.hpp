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

#ifndef CLA3P_LRA_XXMATRIX_HPP_
#define CLA3P_LRA_XXMATRIX_HPP_

/**
 * @file
 */

#include <string>
#include <iostream>

#include "cla3p/types/integer.hpp"
#include "cla3p/generic/meta2d.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
namespace lra {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense real matrix class.
 * @details Represents a low-rank approximation matrix in factored form @f$ A \cdot B^* @f$.
 * @tparam T_Matrix The underlying matrix type (e.g., dns::XxMatrix<T_Scalar>).
 */
template <typename T_Matrix>
class XxMatrix : public Meta2D<int_t> {

	private:
		using T_Scalar = typename T_Matrix::value_type;

	public:

		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Creates an empty low-rank matrix with no allocated memory.
		 */
		XxMatrix();

		/**
		 * @brief Dimension constructor.
		 * @details Constructs a (nr x nc) low-rank matrix with rank @p k and uninitialized values.
		 * @param[in] nr The number of matrix rows.
		 * @param[in] nc The number of matrix columns.
		 * @param[in] k The low matrix rank.
		 */
		explicit XxMatrix(int_t nr, int_t nc, int_t k);

		/**
		 * @brief Copy constructor.
		 * @details Creates a new low-rank matrix by copying another low-rank matrix.
		 * @param[in] other The low-rank matrix to copy.
		 */
		XxMatrix(const XxMatrix<T_Matrix>& other);

		/**
		 * @brief Move constructor.
		 * @details Creates a new low-rank matrix by moving resources from another low-rank matrix.
		 * @param[in] other The low-rank matrix to move from.
		 */
		XxMatrix(XxMatrix<T_Matrix>&& other);

		/**
		 * @brief Destructor.
		 * @details Destroys the low-rank matrix and releases allocated memory.
		 */
		~XxMatrix();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another low-rank matrix to this low-rank matrix.
		 * @param[in] other The low-rank matrix to copy.
		 * @return Reference to this low-rank matrix.
		 */
		XxMatrix<T_Matrix>& operator=(const XxMatrix<T_Matrix>& other);

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another low-rank matrix to this low-rank matrix.
		 * @param[in] other The low-rank matrix to move from.
		 * @return Reference to this low-rank matrix.
		 */
		XxMatrix<T_Matrix>& operator=(XxMatrix<T_Matrix>&& other);

		/**
		 * @brief Fill operator.
		 * @details Fills all elements in both matrices A and B with the specified value.
		 * @param[in] val The scalar value to fill with.
		 */
		void operator=(T_Scalar val);

		/** @} */

		/**
		 * @name Arguments
		 * @{
		 */

		/**
		 * @brief The matrix rows.
		 * @details Returns the number of rows in the low-rank matrix.
		 * @return The number of rows in @c (*this) i.e. number of rows in @c (*this).A().
		 */
		int_t nrows() const override;

		/**
		 * @brief The matrix columns.
		 * @details Returns the number of columns in the low-rank matrix.
		 * @return The number of columns in @c (*this) i.e. number of rows in @c (*this).B().
		 */
		int_t ncols() const override;

		/**
		 * @brief The matrix rank.
		 * @details Returns the rank of the low-rank approximation.
		 * @return The matrix rank i.e. number of columns in @c (*this).A() and @c (*this).B().
		 */
		int_t nrank() const;

		/**
		 * @copydoc cla3p::lra::XxMatrix::A()
		 */
		const T_Matrix& A() const;

		/**
		 * @brief The left matrix of the low-rank product.
		 * @details Returns a reference to the matrix @p A of the product @f$ A \cdot B^* @f$ that is the low-rank approximation.
		 * @return Reference to the left matrix.
		 */
		T_Matrix& A();

		/**
		 * @copydoc cla3p::lra::XxMatrix::B()
		 */
		const T_Matrix& B() const;

		/**
		 * @brief The right matrix of the low-rank product.
		 * @details Returns a reference to the matrix @p B of the product @f$ A \cdot B^* @f$ that is the low-rank approximation.
		 * @return Reference to the right matrix.
		 */
		T_Matrix& B();

		/** @} */

		/**
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Clear the low-rank matrix.
		 * @details Releases all memory and resets the matrix to an empty state.
		 */
		void clear();

		/**
		 * @brief Fill all elements with a value.
		 * @details Sets all elements in both matrices A and B to the specified scalar value.
		 * @param[in] val The scalar value to fill with.
		 */
		void fill(T_Scalar val);

		/**
		 * @brief Scale the low-rank matrix in-place.
		 * @details Multiplies all elements of matrix A by a scalar value.
		 * @param[in] val The scalar value to multiply by.
		 */
		void iscale(T_Scalar val);

		/**
		 * @brief Get information about the low-rank matrix.
		 * @details Returns a string containing information about the matrix's dimensions and rank.
		 * @param[in] header Optional header string to prepend to the information.
		 * @return A string containing matrix information.
		 */
		std::string info(const std::string& header = "") const;

		/** @} */

	private:
		T_Matrix m_A;
		T_Matrix m_B;

		void moveFrom(XxMatrix<T_Matrix>& other);
		void copyFromExisting(const XxMatrix<T_Matrix>& other);
};

/*-------------------------------------------------*/
} // namespace lra
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 * @details Outputs both the left matrix @p A and right matrix @p B of the low-rank factorization.
 * @tparam T_Matrix The underlying matrix type (e.g., dns::XxMatrix<T_Scalar>).
 * @param[in,out] os The output stream.
 * @param[in] mat The low-rank matrix to output.
 * @return The output stream.
 */
template <typename T_Matrix>
std::ostream& operator<<(std::ostream& os, const cla3p::lra::XxMatrix<T_Matrix>& mat)
{
	os << "A:\n";
	os << mat.A();
	os << "B:\n";
	os << mat.B();
	return os;
}

#endif // CLA3P_LRA_XXMATRIX_HPP_
