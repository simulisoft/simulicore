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

#ifndef CULITE_DNS_XXMATRIX_HPP_
#define CULITE_DNS_XXMATRIX_HPP_

/**
 * @file
 */

//#include <string>

#include <cla3p/generic/guard.hpp>
#include <cla3p/generic/matrix_meta.hpp>
#include <cla3p/dense/dns_xxmatrix.hpp>

#include "culite/dense/dns_xxcontainer.hpp"
#include "culite/dense/dns_xxvector.hpp"

#include "culite/virtuals/virtual_transpose.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The device dense matrix class.
 * @details Represents a dense matrix stored on the GPU device.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 */
template <typename T_Scalar>
class XxMatrix : public ::cla3p::MatrixMeta<int_t>, public XxContainer<T_Scalar> {

	private:
		using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
		using T_Cla3pScalar = typename TypeTraits<T_Scalar>::cla3p_type;

	public:
	
		/**
		 * @name Constructors
		 * @{
		 */

		/**
		 * @brief Default constructor.
		 * @details Creates an empty device matrix with no allocated device memory.
		 */
		XxMatrix();

		/**
		 * @brief Dimension constructor.
		 * @details Creates a device matrix of the specified dimensions and allocates device memory.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] pr The matrix property (default: General).
		 */
		explicit XxMatrix(int_t nr, int_t nc, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @brief Auxiliary constructor.
		 * @details Creates a device matrix using existing device memory.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] vals Pointer to existing device memory.
		 * @param[in] ldv The leading dimension.
		 * @param[in] bind If true, the matrix does not take ownership of the memory.
		 * @param[in] pr The matrix property (default: General).
		 */
		explicit XxMatrix(int_t nr, int_t nc, T_Scalar *vals, int_t ldv, bool bind, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/**
		 * @brief Copy constructor.
		 * @details Creates a new device matrix by copying another device matrix.
		 * @param[in] other The device matrix to copy.
		 */
		XxMatrix(const XxMatrix<T_Scalar>& other);

		/**
		 * @brief Move constructor.
		 * @details Creates a new device matrix by moving resources from another device matrix.
		 * @param[in] other The device matrix to move from.
		 */
		XxMatrix(XxMatrix<T_Scalar>&& other);

		/**
		 * @brief Destructor.
		 * @details Destroys the device matrix and releases allocated device memory.
		 */
		~XxMatrix();

		/** @} */

		/**
		 * @name Operators
		 * @{
		 */

		/**
		 * @brief Copy assignment operator.
		 * @details Copies the contents of another device matrix to this device matrix.
		 * @param[in] other The device matrix to copy.
		 * @return Reference to this device matrix.
		 */
		XxMatrix<T_Scalar>& operator=(const XxMatrix<T_Scalar>& other);

		/**
		 * @brief Move assignment operator.
		 * @details Moves resources from another device matrix to this device matrix.
		 * @param[in] other The device matrix to move from.
		 * @return Reference to this device matrix.
		 */
		XxMatrix<T_Scalar>& operator=(XxMatrix<T_Scalar>&& other);

		/**
		 * @brief Unary negation operator.
		 * @details Returns a negated copy of the device matrix.
		 * @return A device matrix containing the negated elements.
		 */
		XxMatrix<T_Scalar> operator-() const; // TODO: use virtuals

		/** @} */

		/** 
		 * @name Arguments
		 * @{
		 */

		/**
		 * @brief Get the leading dimension.
		 * @details Returns the leading dimension of the device matrix (distance between consecutive columns).
		 * @return The leading dimension.
		 */
		int_t ld() const;

		/** @} */

		/**
		 * @name Public Member Functions
		 * @{
		 */

		/**
		 * @brief Clear the device matrix.
		 * @details Releases all device memory and resets the matrix to an empty state.
		 */
		void clear();

		/**
		 * @brief Create a deep copy of the device matrix.
		 * @details Creates a new device matrix with its own device memory allocation
		 *          and copies all elements from this matrix to the new matrix.
		 * @return A new device matrix containing a copy of this matrix's data.
		 */
		XxMatrix<T_Scalar> copy() const;

		/**
		 * @brief Create a reference copy (shallow copy) of the device matrix.
		 * @details Creates a new device matrix object that references the same device memory
		 *          as this matrix. Changes to either matrix will affect both.
		 * @return A device matrix that shares device memory with this matrix.
		 */
		XxMatrix<T_Scalar> rcopy();

		/**
		 * @brief Create a guarded reference copy (shallow copy) of the device matrix.
		 * @details Creates a guarded device matrix object that references the same device memory
		 *          as this matrix. The guard ensures the reference is read-only.
		 * @return A guarded device matrix that shares device memory with this matrix.
		 */
		::cla3p::Guard<XxMatrix<T_Scalar>> rcopy() const;

		/**
		 * @brief Move the device matrix's resources.
		 * @details Transfers ownership of the device memory to a new matrix object,
		 *          leaving this matrix in an empty state.
		 * @return A new device matrix containing this matrix's device memory.
		 */
		XxMatrix<T_Scalar> move();

		/**
		 * @brief Get information about the device matrix.
		 * @details Returns a string containing information about the matrix's dimensions and properties.
		 * @param[in] header Optional header string to prepend to the information.
		 * @return A string containing matrix information.
		 */
		std::string info(const std::string& header = "") const;

		/**
		 * @brief Scale the device matrix in-place.
		 * @details Multiplies all elements of the device matrix by a scalar value.
		 * @param[in] val The scalar value to multiply by.
		 */
		void iscale(T_Scalar val);

		/**
		 * @brief Transpose the device matrix.
		 * @details Returns a transposed copy of the device matrix.
		 * @return A device matrix containing the transposed elements.
		 */
		alias::VirtualTrans_dns<T_Scalar> transpose() const;

		/**
		 * @brief Conjugate transpose the device matrix.
		 * @details Returns a conjugate transposed copy of the device matrix.
		 * @return A device matrix containing the conjugate transposed elements.
		 */
		alias::VirtualTrans_dns<T_Scalar> ctranspose() const;

		/**
		 * @brief Compute the complex conjugate.
		 * @details Returns a device matrix containing the complex conjugate of each element.
		 * @return A device matrix with conjugated elements.
		 */
		XxMatrix<T_Scalar> conjugate() const; // TODO: use virtuals

		/**
		 * @brief Conjugate the device matrix in-place.
		 * @details Replaces all elements with their complex conjugates.
		 */
		void iconjugate();

		/**
		 * @brief Compute the 1-norm.
		 * @details Computes the maximum absolute column sum.
		 * @return The 1-norm of the device matrix.
		 */
		T_RScalar normOne() const;

		/**
		 * @brief Compute the infinity norm.
		 * @details Computes the maximum absolute row sum.
		 * @return The infinity norm of the device matrix.
		 */
		T_RScalar normInf() const;

		/**
		 * @brief Compute the maximum norm.
		 * @details Computes the maximum absolute value of all elements.
		 * @return The maximum norm of the device matrix.
		 */
		T_RScalar normMax() const;

		/**
		 * @brief Compute the Frobenius norm.
		 * @details Computes the square root of the sum of squared absolute values of all elements.
		 * @return The Frobenius norm of the device matrix.
		 */
		T_RScalar normFro() const;

		/**
		 * @brief Extract a block as a new device matrix.
		 * @details Creates a new device matrix containing a deep copy of a contiguous block of elements.
		 * @param[in] ibgn The starting row index of the block.
		 * @param[in] jbgn The starting column index of the block.
		 * @param[in] ni The number of rows in the block.
		 * @param[in] nj The number of columns in the block.
		 * @return A new device matrix containing the block's data.
		 */
		XxMatrix<T_Scalar> block(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/**
		 * @brief Extract a reference block.
		 * @details Creates a device matrix that references a contiguous block of this matrix's device memory.
		 * @param[in] ibgn The starting row index of the block.
		 * @param[in] jbgn The starting column index of the block.
		 * @param[in] ni The number of rows in the block.
		 * @param[in] nj The number of columns in the block.
		 * @return A device matrix that shares device memory with this matrix.
		 */
		XxMatrix<T_Scalar> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj);

		/**
		 * @brief Extract a guarded reference block.
		 * @details Creates a guarded device matrix that references a contiguous block of this matrix's device memory.
		 * @param[in] ibgn The starting row index of the block.
		 * @param[in] jbgn The starting column index of the block.
		 * @param[in] ni The number of rows in the block.
		 * @param[in] nj The number of columns in the block.
		 * @return A guarded device matrix that shares device memory with this matrix.
		 */
		::cla3p::Guard<XxMatrix<T_Scalar>> rblock(int_t ibgn, int_t jbgn, int_t ni, int_t nj) const;

		/**
		 * @brief Set a block of elements.
		 * @details Copies elements from the source device matrix to a contiguous block of this matrix.
		 * @param[in] ibgn The starting row index where the block will be written.
		 * @param[in] jbgn The starting column index where the block will be written.
		 * @param[in] src The source device matrix to copy from.
		 */
		void setBlock(int_t ibgn, int_t jbgn, const XxMatrix<T_Scalar>& src);

		/**
		 * @brief Extract a column as a new device vector.
		 * @details Creates a new device vector containing a deep copy of a column.
		 * @param[in] j The column index.
		 * @return A new device vector containing the column's data.
		 */
		XxVector<T_Scalar> column(int_t j) const;

		/**
		 * @brief Extract a column as a reference.
		 * @details Creates a device vector that references a column of this matrix's device memory.
		 * @param[in] j The column index.
		 * @return A device vector that shares device memory with this matrix.
		 */
		XxVector<T_Scalar> rcolumn(int_t j);

		/**
		 * @brief Extract a column as a guarded reference.
		 * @details Creates a guarded device vector that references a column of this matrix's device memory.
		 * @param[in] j The column index.
		 * @return A guarded device vector that shares device memory with this matrix.
		 */
		::cla3p::Guard<XxVector<T_Scalar>> rcolumn(int_t j) const;

		/**
		 * @brief Extract a row as a new device matrix.
		 * @details Creates a new device matrix containing a deep copy of a row.
		 * @param[in] i The row index.
		 * @return A new device matrix containing the row's data.
		 */
		XxMatrix<T_Scalar> row(int_t i) const;

		/**
		 * @brief Extract a row as a reference.
		 * @details Creates a device matrix that references a row of this matrix's device memory.
		 * @param[in] i The row index.
		 * @return A device matrix that shares device memory with this matrix.
		 */
		XxMatrix<T_Scalar> rrow(int_t i);

		/**
		 * @brief Extract a row as a guarded reference.
		 * @details Creates a guarded device matrix that references a row of this matrix's device memory.
		 * @param[in] i The row index.
		 * @return A guarded device matrix that shares device memory with this matrix.
		 */
		::cla3p::Guard<XxMatrix<T_Scalar>> rrow(int_t i) const;

		/* TODO: use virtuals
		 * @copydoc standard_matrix_docs::rrowvec()
		 */
		//VirtualRowvec<T_Scalar> rrowvec(int_t i) const;

		/**
		 * @brief Copies the device matrix to a host matrix.
		 * @details Performs a device-to-host copy operation.
		 *          If @p dest is empty, it is resized to match the dimensions and properties of the device matrix.
		 *          Otherwise, @p dest must be compatible with the device matrix.
		 * @param[out] dest The destination host matrix.
		 */
		void copyToHost(::cla3p::dns::XxMatrix<T_Cla3pScalar>& dest) const;

		/**
		 * @brief Copies a host matrix to the device matrix.
		 * @details Performs a host-to-device copy operation.
		 *          If the device matrix is empty, it is resized to match the dimensions and properties of @p src.
		 *          Otherwise, the device matrix must be compatible with @p src.
		 * @param[in] src The source host matrix.
		 */
		void copyFromHost(const ::cla3p::dns::XxMatrix<T_Cla3pScalar>& src);

		/** @} */

		/** 
		 * @name Creators/Generators
		 * @{
		 */

		/**
		 * @brief Create a view of existing device memory.
		 * @details Creates a guarded device matrix that references existing device memory
		 *          without taking ownership. The memory must remain valid for the lifetime
		 *          of the returned view.
		 * @param[in] nr The number of rows.
		 * @param[in] nc The number of columns.
		 * @param[in] vals Pointer to the device memory.
		 * @param[in] ldv The leading dimension.
		 * @param[in] pr The matrix property (default: General).
		 * @return A guarded device matrix that views the specified device memory.
		 */
		static ::cla3p::Guard<XxMatrix<T_Scalar>> view(int_t nr, int_t nc, const T_Scalar *vals, int_t ldv, const ::cla3p::Property& pr = ::cla3p::Property::General());

		/** @} */

	private:
		int_t m_ld;

		void defaults();

		void setLd(int_t ld);

		XxMatrix<T_Scalar>& moveFrom(XxMatrix<T_Scalar>& other);
		XxMatrix<T_Scalar>& copyFromExisting(const XxMatrix<T_Scalar>& other);
		void checker() const;
};

/*-------------------------------------------------*/
} // namespace dns
} // namespace culite
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_device2host_operators
 * @brief Copies the contents of src to dest.
 * @details Performs a device-to-host copy from @p src to @p dest.
 *          If @p dest is empty, it is resized to match @p src, otherwise it must be compatible with @p src.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] src The source device matrix.
 * @param[out] dest The destination host matrix.
 */
template <typename T_Scalar>
void operator>>(const culite::dns::XxMatrix<T_Scalar>& src,
	            ::cla3p::dns::XxMatrix<typename culite::TypeTraits<T_Scalar>::cla3p_type>& dest)
{
	src.copyToHost(dest);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_host2device_operators
 * @brief Copies the contents of src to dest.
 * @details Performs a host-to-device copy from @p src to @p dest.
 *          If @p dest is empty, it is resized to match @p src, otherwise it must be compatible with @p src.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in] src The source host matrix.
 * @param[out] dest The destination device matrix.
 */
template <typename T_Scalar>
void operator>>(const ::cla3p::dns::XxMatrix<typename culite::TypeTraits<T_Scalar>::cla3p_type>& src,
	            culite::dns::XxMatrix<T_Scalar>& dest)
{
	dest.copyFromHost(src);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 * @details Copies the device matrix to host and writes it to the output stream.
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 * @param[in,out] os The output stream.
 * @param[in] mat The device matrix to output.
 * @return The output stream.
 */
template <typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const culite::dns::XxMatrix<T_Scalar>& mat)
{
	::cla3p::dns::XxMatrix<typename culite::TypeTraits<T_Scalar>::cla3p_type> hostMat;
	mat >> hostMat;
	os << hostMat;
	return os;
}

/*-------------------------------------------------*/

#endif // CULITE_DNS_XXMATRIX_HPP_
