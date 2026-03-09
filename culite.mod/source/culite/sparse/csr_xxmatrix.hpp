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

#ifndef CULITE_CSR_XXMATRIX_HPP_
#define CULITE_CSR_XXMATRIX_HPP_

/**
 * @file
 */

#include <ostream>
#include <string>

#include <cla3p/sparse/csr_xxmatrix.hpp>

#include "culite/sparse/csx_xxcontainer.hpp"

#include "culite/virtuals/virtual_expression.hpp"
#include "culite/virtuals/virtual_object.hpp"
#include "culite/virtuals/virtual_transpose.hpp"
#include "culite/virtuals/virtual_conjugate.hpp"
#include "culite/virtuals/virtual_scale.hpp"

// forwards
#include "culite/types/cla3p_forwards.hpp"
#include "culite/generic/cla3p_forwards.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace csr {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The device sparse matrix class (compressed sparse row format).
 * @details Represents a sparse matrix stored on the GPU device in CSR format.
 * @tparam T_Int The integer type for indexing (e.g., int32_t, int64_t).
 * @tparam T_Scalar The scalar type (e.g., float, double, complex).
 */
template <typename T_Int, typename T_Scalar>
class XxMatrix : public MatrixMeta<T_Int>, public csx::XxContainer<T_Int,T_Scalar> {

    private:
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
        using T_Cla3pInt = typename TypeTraits<T_Int>::cla3p_type;
        using T_Cla3pScalar = typename TypeTraits<T_Scalar>::cla3p_type;

    public:

        /**
         * @name Virtual Convertors
         * @{
         */

        template <typename T_Virtual>
        XxMatrix(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v) { evaluateFrom(v); }

        template <typename T_Virtual>
        XxMatrix<T_Int,T_Scalar>& operator=(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v) { return evaluateFrom(v); }

        alias::VirtualObj_csr<T_Int,T_Scalar> virtualize() const { return alias::VirtualObj_csr<T_Int,T_Scalar>(*this); }

        /** @} */

        /**
         * @name Constructors
         * @{
         */

        /**
         * @brief Default constructor.
         * @details Creates an empty device sparse matrix with no allocated device memory.
         */
        XxMatrix();

        /**
         * @brief Dimension constructor.
         * @details Creates a device sparse matrix of the specified dimensions in CSR format and allocates device memory.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] nz The number of non-zero elements.
         * @param[in] pr The matrix property (default: General).
         */
        explicit XxMatrix(T_Int nr, T_Int nc, T_Int nz, const Property& pr = Property::General());

        /**
         * @brief Auxiliary constructor.
         * @details Creates a device sparse matrix using existing device memory CSR format arrays.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] rptr Pointer to existing device row pointer array.
         * @param[in] cidx Pointer to existing device column index array.
         * @param[in] vals Pointer to existing device values array.
         * @param[in] bind If true, the matrix takes ownership of the device memory.
         * @param[in] pr The matrix property (default: General).
         */
        explicit XxMatrix(T_Int nr, T_Int nc, T_Int *rptr, T_Int *cidx, T_Scalar *vals, bool bind, 
                          const Property& pr = Property::General());

        /**
         * @brief Copy constructor.
         * @details Creates a new device sparse matrix by copying another device sparse matrix.
         * @param[in] other The device sparse matrix to copy.
         */
        XxMatrix(const XxMatrix<T_Int,T_Scalar>& other);

        /**
         * @brief Move constructor.
         * @details Creates a new device sparse matrix by moving resources from another device sparse matrix.
         * @param[in] other The device sparse matrix to move from.
         */
        XxMatrix(XxMatrix<T_Int,T_Scalar>&& other);

        /**
         * @brief Destructor.
         * @details Destroys the device sparse matrix and releases allocated device memory.
         */
        ~XxMatrix();

        /** @} */

        /**
         * @name Operators
         * @{
         */

        /**
         * @brief Copy assignment operator.
         * @details Copies the contents of another device sparse matrix to this device sparse matrix.
         * @param[in] other The device sparse matrix to copy.
         * @return Reference to this device sparse matrix.
         */
        XxMatrix<T_Int,T_Scalar>& operator=(const XxMatrix<T_Int,T_Scalar>& other);

        /**
         * @brief Move assignment operator.
         * @details Moves resources from another device sparse matrix to this device sparse matrix.
         * @param[in] other The device sparse matrix to move from.
         * @return Reference to this device sparse matrix.
         */
        XxMatrix<T_Int,T_Scalar>& operator=(XxMatrix<T_Int,T_Scalar>&& other);

        /**
         * @brief Unary negation operator.
         * @details Returns a negated copy of the device sparse matrix.
         * @return A device sparse matrix containing the negated elements.
         */
        alias::VirtualScal_csr<T_Int,T_Scalar> operator-() const;

        /** @} */

        /** 
         * @name Arguments
         * @{
         */

        /**
         * @brief Access the row pointer array.
         * @details Returns a pointer to the row pointer array (CSR format).
         * @return Pointer to the row pointer array.
         */
        T_Int* rowptr();

        /**
         * @brief Access the row pointer array.
         * @details Returns a pointer to the row pointer array (CSR format).
         * @return Pointer to the row pointer array.
         */
        const T_Int* rowptr() const;

        /**
         * @brief Access the column index array.
         * @details Returns a pointer to the column index array (CSR format).
         * @return Pointer to the column index array.
         */
        T_Int* colidx();

        /**
         * @brief Access the column index array.
         * @details Returns a pointer to the column index array (CSR format).
         * @return Pointer to the column index array.
         */
        const T_Int* colidx() const;

        /**
         * @brief Get the number of non-zero elements.
         * @details Returns the number of stored non-zero elements in the device sparse matrix.
         * @return The number of non-zero elements.
         */
        T_Int nnz() const;

        /** @} */

        /** 
         * @name Public Member Functions
         * @{
         */

        /**
         * @brief Clear the device sparse matrix.
         * @details Releases all device memory and resets the matrix to an empty state.
         */
        void clear();

        /**
         * @brief Get information about the sparse matrix.
         * @details Returns a string containing information about the matrix's dimensions, non-zeros, and properties.
         * @param[in] header Optional header string to prepend to the information.
         * @return A string containing matrix information.
         */
        std::string info(const std::string& header = "") const;

        /**
         * @brief Create a deep copy of the device sparse matrix.
         * @details Creates a new device sparse matrix with its own device memory allocation
         *          and copies all elements from this matrix to the new matrix.
         * @return A new device sparse matrix containing a copy of this matrix's data.
         */
        XxMatrix<T_Int,T_Scalar> copy() const;

        /**
         * @brief Create a reference copy (shallow copy) of the device sparse matrix.
         * @details Creates a new device sparse matrix object that references the same device memory
         *          as this matrix. Changes to either matrix will affect both.
         * @return A device sparse matrix that shares device memory with this matrix.
         */
        XxMatrix<T_Int,T_Scalar> rcopy();

        /**
         * @brief Create a guarded reference copy (shallow copy) of the device sparse matrix.
         * @details Creates a guarded device sparse matrix object that references the same device memory
         *          as this matrix. The guard ensures the reference is read-only.
         * @return A guarded device sparse matrix that shares device memory with this matrix.
         */
        Guard<XxMatrix<T_Int,T_Scalar>> rcopy() const;

        /**
         * @brief Move the device sparse matrix's resources.
         * @details Transfers ownership of the device memory to a new matrix object,
         *          leaving this matrix in an empty state.
         * @return A new device sparse matrix containing this matrix's device memory.
         */
        XxMatrix<T_Int,T_Scalar> move();

        /**
         * @brief Scale the device sparse matrix in-place.
         * @details Multiplies all non-zero elements of the device sparse matrix by a scalar value.
         * @param[in] val The scalar value to multiply by.
         */
        void iscale(T_Scalar val);

        /**
         * @brief Compute the transpose.
         * @details Returns a virtual expression representing the transpose of the device sparse matrix.
         *          Rows and columns are swapped.
         * @return A virtual expression for the transposed matrix.
         * @warning Explicit calculation is not supported. Evaluating the returned expression on an existing matrix will throw an exception.
         */
        alias::VirtualTrans_csr<T_Int,T_Scalar> transpose() const;

        /**
         * @brief Compute the conjugate transpose (Hermitian transpose).
         * @details Returns a virtual expression representing the conjugate transpose of the device sparse matrix.
         *          Rows and columns are swapped, and complex elements are conjugated.
         *          For real matrices, this is equivalent to transpose().
         * @return A virtual expression for the conjugate transposed matrix.
         * @warning Explicit calculation is not supported. Evaluating the returned expression on an existing matrix will throw an exception.
         */
        alias::VirtualTrans_csr<T_Int,T_Scalar> ctranspose() const;

        /**
         * @brief Compute the complex conjugate.
         * @details Returns a device matrix containing the complex conjugate of each non-zero element.
         * @return A device sparse matrix with conjugated elements.
         */
        alias::VirtualConj_csr<T_Int,T_Scalar> conjugate() const;

        /**
         * @brief Conjugate the device sparse matrix in-place.
         * @details Replaces all non-zero elements with their complex conjugates.
         */
        void iconjugate();

        /**
         * @brief Copy the device sparse matrix to host memory.
         * @details Transfers the device sparse matrix data to a host sparse matrix.
         * @param[out] dest The host sparse matrix destination.
         */
        void copyToHost(::cla3p::csr::XxMatrix<T_Cla3pInt, T_Cla3pScalar>& dest) const;

        /** 
         * @brief Copy a host sparse matrix to this device sparse matrix.
         * @details Transfers data from a host sparse matrix to this device sparse matrix.
         * @param[in] src The host sparse matrix source.
         */
        void copyFromHost(const ::cla3p::csr::XxMatrix<T_Cla3pInt, T_Cla3pScalar>& src);

        /** @} */

        /** 
         * @name Creators/Generators
         * @{
         */

        /**
         * @brief Create a view of existing device CSR arrays.
         * @details Creates a guarded device sparse matrix that references existing device CSR format memory
         *          without taking ownership. The device memory must remain valid for the lifetime
         *          of the returned view.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] rptr Pointer to the device row pointer array.
         * @param[in] cidx Pointer to the device column index array.
         * @param[in] vals Pointer to the device values array.
         * @param[in] pr The matrix property (default: General).
         * @return A guarded device sparse matrix that views the specified device CSR arrays.
         */
        static Guard<XxMatrix<T_Int,T_Scalar>> 
        view(T_Int nr, T_Int nc, 
             const T_Int *rptr, const T_Int *cidx, const T_Scalar *vals, 
             const Property& pr = Property::General());

        /** @} */

    private:
        XxMatrix<T_Int,T_Scalar>& copyFromExisting(const XxMatrix<T_Int,T_Scalar>& other);
        XxMatrix<T_Int,T_Scalar>& moveFrom(XxMatrix<T_Int,T_Scalar>& other);
        void checker() const;

        template <typename T_Virtual>
        XxMatrix<T_Int,T_Scalar>& evaluateFrom(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v)
        {
            if(*this) {
                v.evaluateOnExisting(*this);
            } else {
                v.evaluateOnNew(*this);
            }
            return *this;
        }
};

/*-------------------------------------------------*/
} // namespace csr
} // namespace culite
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Stream operator for copying device sparse matrix to host.
 * @details Transfers the contents of a device sparse matrix to a host sparse matrix.
 *          This operator enables convenient syntax for device-to-host memory transfers.
 * @tparam T_Int The integer type for indexing.
 * @tparam T_Scalar The scalar type.
 * @param[in] src The source device sparse matrix.
 * @param[out] dest The destination host sparse matrix.
 */
template <typename T_Int, typename T_Scalar>
void operator>>(const culite::csr::XxMatrix<T_Int,T_Scalar>& src,
                cla3p::csr::XxMatrix<typename culite::TypeTraits<T_Int>::cla3p_type,
                                     typename culite::TypeTraits<T_Scalar>::cla3p_type>& dest)
{
    src.copyToHost(dest);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
* @brief Stream operator for copying host sparse matrix to device.
* @details Transfers the contents of a host sparse matrix to a device sparse matrix.
*          This operator enables convenient syntax for host-to-device memory transfers.
* @tparam T_Int The integer type for indexing.
* @tparam T_Scalar The scalar type.
* @param[in] src The source host sparse matrix.
* @param[out] dest The destination device sparse matrix.
 */
template <typename T_Int, typename T_Scalar>
void operator>>(const cla3p::csr::XxMatrix<typename culite::TypeTraits<T_Int>::cla3p_type,
                                           typename culite::TypeTraits<T_Scalar>::cla3p_type>& src,
                culite::csr::XxMatrix<T_Int,T_Scalar>& dest)
{
    dest.copyFromHost(src);
}

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 */
template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const culite::csr::XxMatrix<T_Int,T_Scalar>& mat)
{
    cla3p::csr::XxMatrix<typename culite::TypeTraits<T_Int>::cla3p_type,
                         typename culite::TypeTraits<T_Scalar>::cla3p_type> hostMat;
    mat >> hostMat;
    os << hostMat;
    return os;
}

/*-------------------------------------------------*/

#endif // CULITE_CSR_XXMATRIX_HPP_