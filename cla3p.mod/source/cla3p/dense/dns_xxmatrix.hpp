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

#ifndef CLA3P_DNS_XXMATRIX_HPP_
#define CLA3P_DNS_XXMATRIX_HPP_

/**
 * @file
 */

#include <string>

#include "cla3p/generic/guard.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/dense/dns_xxcontainer.hpp"
#include "cla3p/dense/dns_xxvector.hpp"

#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_rowvec.hpp"
#include "cla3p/virtuals/virtual_transpose.hpp"
#include "cla3p/virtuals/virtual_conjugate.hpp"
#include "cla3p/virtuals/virtual_scale.hpp"

/*-------------------------------------------------*/

namespace cla3p { namespace prm { template <typename T_Int> class PxMatrix; } }

/*-------------------------------------------------*/
namespace cla3p { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The dense matrix class.
 * @details Represents a dense matrix stored in column-major format.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 */
template <typename T_Scalar>
class XxMatrix : public MatrixMeta<int_t>, public XxContainer<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

	public:

		//
		// Convertors
		// Move convertors intentionally left as non-explicit
		//
		template <typename T_Virtual>
		XxMatrix(const alias::VirtualExpr_dns<T_Scalar,T_Virtual>& v) { operator=(v); }
		template <typename T_Virtual>
		XxMatrix<T_Scalar>& operator=(const alias::VirtualExpr_dns<T_Scalar,T_Virtual>& v) { return evaluateFrom(v); }

		//template <typename T_Virtual>
		//XxMatrix(alias::VirtualExpr_dns<T_Scalar,T_Virtual>&& v) { operator=(std::move(v)); }
		//template <typename T_Virtual>
		//XxMatrix<T_Scalar>& operator=(alias::VirtualExpr_dns<T_Scalar,T_Virtual>&& v) { return evaluateFrom(v); }

		XxMatrix(const VirtualRowvec<T_Scalar>& rv) { operator=(rv); }
        XxMatrix<T_Scalar>& operator=(const VirtualRowvec<T_Scalar>& rv) { return evaluateFrom(rv); }
        
		//XxMatrix(VirtualRowvec<T_Scalar>&& rv) { operator=(std::move(rv)); }
		//XxMatrix<T_Scalar>& operator=(VirtualRowvec<T_Scalar>&& rv) { return evaluateFrom(rv);  }
		
		alias::VirtualObj_dns<T_Scalar> virtualize() const { return alias::VirtualObj_dns<T_Scalar>(*this); }
	
		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Creates an empty matrix with no allocated memory.
		 */
		XxMatrix();

		/**
		 * @brief Dimension constructor.
		 * @details Creates a matrix of the specified dimensions and allocates memory.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] pr The matrix property (default: General).
		 */
		explicit XxMatrix(int_t nr, int_t nc, const Property& pr = Property::General());

		/**
		 * @brief Auxiliary constructor.
		 * @details Creates a matrix using existing memory.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] vals Pointer to existing memory.
		 * @param[in] ldv The leading dimension.
		 * @param[in] bind If true, the matrix does not take ownership of the memory.
		 * @param[in] pr The matrix property (default: General).
		 */
		explicit XxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const Property& pr = Property::General());

		/**
		 * @brief Copy constructor.
		 * @details Creates a new matrix by copying another matrix.
		 * @param[in] other The matrix to copy.
		 */
		XxMatrix(const XxMatrix<T_Scalar>& other);

		/**
		 * @brief Move constructor.
		 * @details Creates a new matrix by moving resources from another matrix.
		 * @param[in] other The matrix to move from.
		 */
		XxMatrix(XxMatrix<T_Scalar>&& other);

		/**
		 * @brief Destructor.
		 * @details Destroys the matrix and releases allocated memory.
		 */
		~XxMatrix();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another matrix to this matrix.
		 * @param[in] other The matrix to copy.
		 * @return Reference to this matrix.
		 */
		XxMatrix<T_Scalar>& operator=(const XxMatrix<T_Scalar>& other);

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another matrix to this matrix.
		 * @param[in] other The matrix to move from.
		 * @return Reference to this matrix.
		 */
		XxMatrix<T_Scalar>& operator=(XxMatrix<T_Scalar>&& other);

		/**
		 * @brief Element access operator.
		 * @details Accesses the element at the specified row and column.
		 * @param[in] i The row index (0-based).
		 * @param[in] j The column index (0-based).
		 * @return Reference to the element at position (i,j).
		 */
		T_Scalar& operator()(int_t i, int_t j);

		/**
		 * @brief Element access operator.
		 * @details Accesses the element at the specified row and column.
		 * @param[in] i The row index (0-based).
		 * @param[in] j The column index (0-based).
		 * @return Const reference to the element at position (i,j).
		 */
		const T_Scalar& operator()(int_t i, int_t j) const;

		/**
		 * @brief Unary negation operator.
		 * @details Returns a negated copy of the matrix.
		 * @return A virtual expression containing the negated elements.
		 */
		alias::VirtualScal_dns<T_Scalar> operator-() const;

		/**
		 * @brief Fill operator.
		 * @details Fills all elements with the specified value.
		 * @param[in] val The scalar value to fill with.
		 */
		void operator=(T_Scalar val);

		/** @} */

		/** 
		 * @name Arguments
		 * @{
		 */

		/**
		 * @brief Get the leading dimension.
		 * @details Returns the leading dimension of the matrix (distance between consecutive columns).
		 * @return The leading dimension.
		 */
		int_t ld() const;

		/** @} */

		/**
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Clear the matrix.
		 * @details Releases all memory and resets the matrix to an empty state.
		 */
		void clear();

		/**
		 * @brief Fill all elements with a value.
		 * @details Sets all elements in the matrix to the specified scalar value.
		 * @param[in] val The scalar value to fill with.
		 */
		void fill(T_Scalar val);

		/**
		 * @brief Create a deep copy of the matrix.
		 * @details Creates a new matrix with its own memory allocation
		 *          and copies all elements from this matrix to the new matrix.
		 * @return A new matrix containing a copy of this matrix's data.
		 */
		XxMatrix<T_Scalar> copy() const;

		/**
		 * @brief Create a reference copy (shallow copy) of the matrix.
		 * @details Creates a new matrix object that references the same memory
		 *          as this matrix. Changes to either matrix will affect both.
		 * @return A matrix that shares memory with this matrix.
		 */
		XxMatrix<T_Scalar> rcopy();

		/**
		 * @brief Create a guarded reference copy (shallow copy) of the matrix.
		 * @details Creates a guarded matrix object that references the same memory
		 *          as this matrix. The guard ensures the reference is read-only.
		 * @return A guarded matrix that shares memory with this matrix.
		 */
		Guard<XxMatrix<T_Scalar>> rcopy() const;

		/**
		 * @brief Move the matrix's resources.
		 * @details Transfers ownership of the memory to a new matrix object,
		 *          leaving this matrix in an empty state.
		 * @return A new matrix containing this matrix's memory.
		 */
		XxMatrix<T_Scalar> move();

		/**
		 * @brief Get information about the matrix.
		 * @details Returns a string containing information about the matrix's dimensions and properties.
		 * @param[in] header Optional header string to prepend to the information.
		 * @return A string containing matrix information.
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @brief Print the matrix to standard output.
		 * @details Prints the matrix elements to standard output with specified precision.
		 * @param[in] prec The output precision (default: 0 for default precision).
		 */
		void print(std::streamsize prec = 0) const;

		/**
		 * @brief Output the matrix to a stream.
		 * @details Writes the matrix elements to the specified output stream with specified precision.
		 * @param[in,out] os The output stream to write to.
		 * @param[in] prec The output precision (default: 0 for default precision).
		 */
		void toStream(std::ostream& os, std::streamsize prec = 0) const;

		/**
		 * @brief Scale the matrix in-place.
		 * @details Multiplies all elements of the matrix by a scalar value.
		 * @param[in] val The scalar value to multiply by.
		 */
		void iscale(T_Scalar val);

		/**
		 * @brief Transpose the matrix.
		 * @details Returns a transposed view of the matrix.
		 * @return A virtual transpose expression.
		 */
		VirtualTranspose<XxMatrix<T_Scalar>> transpose() const;

		/**
		 * @brief Conjugate transpose the matrix.
		 * @details Returns a conjugate transposed view of the matrix.
		 * @return A virtual conjugate transpose expression.
		 */
		VirtualTranspose<XxMatrix<T_Scalar>> ctranspose() const;

		/**
		 * @brief Compute the complex conjugate.
		 * @details Returns a virtual expression containing the complex conjugate of each element.
		 * @return A virtual expression with conjugated elements.
		 */
		alias::VirtualConj_dns<T_Scalar> conjugate() const;

		/**
		 * @brief Conjugate the matrix in-place.
		 * @details Replaces all elements with their complex conjugates.
		 */
		void iconjugate();

		/**
		 * @brief Compute the 1-norm.
		 * @details Computes the maximum absolute column sum.
		 * @return The 1-norm of the matrix.
		 */
		T_RScalar normOne() const;

		/**
		 * @brief Compute the infinity norm.
		 * @details Computes the maximum absolute row sum.
		 * @return The infinity norm of the matrix.
		 */
		T_RScalar normInf() const;

		/**
		 * @brief Compute the maximum norm.
		 * @details Computes the maximum absolute value of all elements.
		 * @return The maximum norm of the matrix.
		 */
		T_RScalar normMax() const;

		/**
		 * @brief Compute the Frobenius norm.
		 * @details Computes the square root of the sum of squared absolute values of all elements.
		 * @return The Frobenius norm of the matrix.
		 */
		T_RScalar normFro() const;

		/**
		 * @brief Convert to general matrix.
		 * @details Creates a new general matrix by copying all relevant elements from a structured matrix.
		 * @return A general matrix containing the full data.
		 */
		XxMatrix<T_Scalar> general() const;

		/**
		 * @brief Convert to general matrix in-place.
		 * @details Converts a structured matrix to general form by filling in all elements.
		 */
		void igeneral();

		/**
		 * @brief Permute rows and columns.
		 * @details Applies left and right permutations to the matrix.
		 * @param[in] P The row permutation matrix.
		 * @param[in] Q The column permutation matrix.
		 * @return A new matrix containing the permuted elements.
		 */
		XxMatrix<T_Scalar> permuteLeftRight(const prm::PxMatrix<int_t>& P, const prm::PxMatrix<int_t>& Q) const;

		/**
		 * @brief Permute rows and columns into a destination.
		 * @details Applies left and right permutations to the matrix and stores the result in the destination matrix.
		 * @param[in] P The row permutation matrix.
		 * @param[in] Q The column permutation matrix.
		 * @param[out] dest The destination matrix for the permuted elements.
		 */
		void permuteLeftRight(const prm::PxMatrix<int_t>& P, const prm::PxMatrix<int_t>& Q, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @brief Permute rows.
		 * @details Applies a left (row) permutation to the matrix.
		 * @param[in] P The row permutation matrix.
		 * @return A new matrix containing the permuted elements.
		 */
		XxMatrix<T_Scalar> permuteLeft(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @brief Permute rows into a destination.
		 * @details Applies a left (row) permutation to the matrix and stores the result in the destination matrix.
		 * @param[in] P The row permutation matrix.
		 * @param[out] dest The destination matrix for the permuted elements.
		 */
		void permuteLeft(const prm::PxMatrix<int_t>& P, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @brief Permute columns.
		 * @details Applies a right (column) permutation to the matrix.
		 * @param[in] Q The column permutation matrix.
		 * @return A new matrix containing the permuted elements.
		 */
		XxMatrix<T_Scalar> permuteRight(const prm::PxMatrix<int_t>& Q) const;

		/**
		 * @brief Permute columns into a destination.
		 * @details Applies a right (column) permutation to the matrix and stores the result in the destination matrix.
		 * @param[in] Q The column permutation matrix.
		 * @param[out] dest The destination matrix for the permuted elements.
		 */
		void permuteRight(const prm::PxMatrix<int_t>& Q, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @brief Permute rows and columns symmetrically.
		 * @details Applies the same permutation to both rows and columns (mirror permutation).
		 * @param[in] P The permutation matrix to apply to both rows and columns.
		 * @return A new matrix containing the permuted elements.
		 */
		XxMatrix<T_Scalar> permuteMirror(const prm::PxMatrix<int_t>& P) const;

		/**
		 * @brief Permute rows and columns symmetrically into a destination.
		 * @details Applies the same permutation to both rows and columns and stores the result in the destination matrix.
		 * @param[in] P The permutation matrix to apply to both rows and columns.
		 * @param[out] dest The destination matrix for the permuted elements.
		 */
		void permuteMirror(const prm::PxMatrix<int_t>& P, XxMatrix<T_Scalar>& dest) const;

		/**
		 * @brief Extract a block as a new matrix.
		 * @details Creates a new matrix containing a deep copy of a contiguous block of elements.
		 * @param[in] ibgn The starting row index of the block.
		 * @param[in] jbgn The starting column index of the block.
		 * @param[in] ni The number of rows in the block.
		 * @param[in] nj The number of columns in the block.
		 * @return A new matrix containing the block's data.
		 */
		XxMatrix<T_Scalar> block(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/**
		 * @brief Extract a reference block.
		 * @details Creates a matrix that references a contiguous block of this matrix's memory.
		 * @param[in] ibgn The starting row index of the block.
		 * @param[in] jbgn The starting column index of the block.
		 * @param[in] ni The number of rows in the block.
		 * @param[in] nj The number of columns in the block.
		 * @return A matrix that shares memory with this matrix.
		 */
		XxMatrix<T_Scalar> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj);

		/**
		 * @brief Extract a guarded reference block.
		 * @details Creates a guarded matrix that references a contiguous block of this matrix's memory.
		 * @param[in] ibgn The starting row index of the block.
		 * @param[in] jbgn The starting column index of the block.
		 * @param[in] ni The number of rows in the block.
		 * @param[in] nj The number of columns in the block.
		 * @return A guarded matrix that shares memory with this matrix.
		 */
		Guard<XxMatrix<T_Scalar>> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/**
		 * @brief Set a block of elements.
		 * @details Copies elements from the source matrix to a contiguous block of this matrix.
		 * @param[in] ibgn The starting row index where the block will be written.
		 * @param[in] jbgn The starting column index where the block will be written.
		 * @param[in] src The source matrix to copy from.
		 */
		void setBlock(int_t ibgn, int_t jbgn, const XxMatrix<T_Scalar>& src);

		/**
		 * @brief Extract a column as a new vector.
		 * @details Creates a new vector containing a deep copy of a column.
		 * @param[in] j The column index.
		 * @return A new vector containing the column's data.
		 */
		XxVector<T_Scalar> column(int_t j) const;

		/**
		 * @brief Extract a column as a reference.
		 * @details Creates a vector that references a column of this matrix's memory.
		 * @param[in] j The column index.
		 * @return A vector that shares memory with this matrix.
		 */
		XxVector<T_Scalar> rcolumn(int_t j);

		/**
		 * @brief Extract a column as a guarded reference.
		 * @details Creates a guarded vector that references a column of this matrix's memory.
		 * @param[in] j The column index.
		 * @return A guarded vector that shares memory with this matrix.
		 */
		Guard<XxVector<T_Scalar>> rcolumn(int_t j) const;

		/**
		 * @brief Extract a row as a new matrix.
		 * @details Creates a new matrix containing a deep copy of a row.
		 * @param[in] i The row index.
		 * @return A new matrix containing the row's data.
		 */
		XxMatrix<T_Scalar> row(int_t i) const;

		/**
		 * @brief Extract a row as a reference.
		 * @details Creates a matrix that references a row of this matrix's memory.
		 * @param[in] i The row index.
		 * @return A matrix that shares memory with this matrix.
		 */
		XxMatrix<T_Scalar> rrow(int_t i);

		/**
		 * @brief Extract a row as a guarded reference.
		 * @details Creates a guarded matrix that references a row of this matrix's memory.
		 * @param[in] i The row index.
		 * @return A guarded matrix that shares memory with this matrix.
		 */
		Guard<XxMatrix<T_Scalar>> rrow(int_t i) const;

		/**
		 * @brief Extract a row as a virtual row vector.
		 * @details Creates a virtual row vector expression that references a row of this matrix's memory.
		 * @param[in] i The row index.
		 * @return A virtual row vector expression.
		 */
		VirtualRowvec<T_Scalar> rrowvec(int_t i) const;

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @brief Create a random matrix.
		 * @details Creates a matrix with random values uniformly distributed in the specified range.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] pr The matrix property (default: General).
		 * @param[in] lo The lower bound of the random values (default: 0).
		 * @param[in] hi The upper bound of the random values (default: 1).
		 * @return A matrix containing random values.
		 */
		static XxMatrix<T_Scalar> random(int_t nr, int_t nc, const Property& pr = Property::General(), 
				T_RScalar lo = T_RScalar(0), 
				T_RScalar hi = T_RScalar(1));

		/**
		 * @brief Create a view of existing memory.
		 * @details Creates a guarded matrix that references existing memory
		 *          without taking ownership. The memory must remain valid for the lifetime
		 *          of the returned view.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] vals Pointer to the memory.
		 * @param[in] ldv The leading dimension.
		 * @param[in] pr The matrix property (default: General).
		 * @return A guarded matrix that views the specified memory.
		 */
		static Guard<XxMatrix<T_Scalar>> view(int_t nr, int_t nc, const T_Scalar *vals, int_t ldv, const Property& pr = Property::General());

		/** @} */

	private:
		int_t m_ld;

		void defaults();

		void setLd(int_t ld);

		XxMatrix<T_Scalar>& moveFrom(XxMatrix<T_Scalar>& other);
		XxMatrix<T_Scalar>& copyFromExisting(const XxMatrix<T_Scalar>& other);
		void checker() const;

	protected:
		template <typename T_Virtual>
		XxMatrix<T_Scalar>& evaluateFrom(const alias::VirtualExpr_dns<T_Scalar,T_Virtual>& v)
		{
			if(*this) {
				v.evaluateOnExisting(*this);
			} else {
				v.evaluateOnNew(*this);
			}
            return *this;
		}

		XxMatrix<T_Scalar>& evaluateFrom(const VirtualRowvec<T_Scalar>& rv)
		{
			if(*this) {
				rv.evaluateOnExisting(*this);
			} else {
				rv.evaluateOnNew(*this);
			}
            return *this;
		}
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 * @details Outputs the matrix to the stream in a human-readable format.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in,out] os The output stream.
 * @param[in] mat The matrix to output.
 * @return The output stream.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::dns::XxMatrix<T_Scalar>& mat)
{
	mat.toStream(os);
	return os;
}

/*-------------------------------------------------*/

#endif // CLA3P_DNS_XXMATRIX_HPP_
