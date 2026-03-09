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

#ifndef CLA3P_CSR_XXMATRIX_HPP_
#define CLA3P_CSR_XXMATRIX_HPP_

/**
 * @file
 */

#include <ostream>
#include <string>

#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/sparse/csx_xxcontainer.hpp"
#include "cla3p/generic/guard.hpp"

#include "cla3p/dense/dns_xxmatrix.hpp"

#include "cla3p/virtuals/virtual_expression.hpp"
#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_transpose.hpp"
#include "cla3p/virtuals/virtual_scale.hpp"
#include "cla3p/virtuals/virtual_strided.hpp"

/*-------------------------------------------------*/

namespace cla3p{ namespace prm { template <typename T_Int> class PxMatrix; } }

/*-------------------------------------------------*/
namespace cla3p { 
namespace csr {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse matrix class (compressed sparse row format).
 */
template <typename T_Int, typename T_Scalar>
class XxMatrix : public MatrixMeta<T_Int>, public csx::XxContainer<T_Int,T_Scalar> {

    private:
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    public:

        /**
         * @name Virtual Convertors
         * @{
         */

        template <typename T_Virtual>
        XxMatrix(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v) { operator=(v); }
        
        template <typename T_Virtual>
        XxMatrix<T_Int,T_Scalar>& operator=(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v) { return evaluateFrom(v); }

        XxMatrix(const alias::VirtualStrided_csr<T_Int,T_Scalar>& v) { evaluateFrom(v); }
        XxMatrix<T_Int,T_Scalar>& operator=(const alias::VirtualStrided_csr<T_Int,T_Scalar>& v) { return evaluateFrom(v); }

        XxMatrix(const alias::GuardedStrided_csr<T_Int,T_Scalar>& gv) { evaluateFrom(gv); }
        XxMatrix<T_Int,T_Scalar>& operator=(const alias::GuardedStrided_csr<T_Int,T_Scalar>& gv) { return evaluateFrom(gv); }

        alias::VirtualObj_csr<T_Int,T_Scalar> virtualize() const { return alias::VirtualObj_csr<T_Int,T_Scalar>(*this); }

        /** @} */

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
         * @details Creates a sparse matrix of the specified dimensions in CSR format and allocates memory.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] nz The number of non-zero elements.
         * @param[in] pr The matrix property (default: General).
         */
        explicit XxMatrix(T_Int nr, T_Int nc, T_Int nz, const Property& pr = Property::General());

        /**
         * @brief Auxiliary constructor.
         * @details Creates a sparse matrix using existing CSR format arrays.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] rptr Pointer to existing row pointer array.
         * @param[in] cidx Pointer to existing column index array.
         * @param[in] vals Pointer to existing values array.
         * @param[in] bind If true, the matrix takes ownership of the memory.
         * @param[in] pr The matrix property (default: General).
         */
        explicit XxMatrix(T_Int nr, T_Int nc, T_Int *rptr, T_Int *cidx, T_Scalar *vals, bool bind, const Property& pr = Property::General());

        /**
         * @brief Copy constructor.
         * @details Creates a new sparse matrix by copying another sparse matrix.
         * @param[in] other The sparse matrix to copy.
         */
        XxMatrix(const XxMatrix<T_Int,T_Scalar>& other);

        /**
         * @brief Move constructor.
         * @details Creates a new sparse matrix by moving resources from another sparse matrix.
         * @param[in] other The sparse matrix to move from.
         */
        XxMatrix(XxMatrix<T_Int,T_Scalar>&& other);

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
        XxMatrix<T_Int,T_Scalar>& operator=(const XxMatrix<T_Int,T_Scalar>& other);

        /**
         * @brief Move assignment operator.
         * @details Moves resources from another sparse matrix to this sparse matrix.
         * @param[in] other The sparse matrix to move from.
         * @return Reference to this sparse matrix.
         */
        XxMatrix<T_Int,T_Scalar>& operator=(XxMatrix<T_Int,T_Scalar>&& other);

        /**
         * @brief Unary negation operator.
         * @details Returns a negated copy of the sparse matrix.
         * @return A virtual expression containing the negated elements.
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
         * @brief Get information about the sparse matrix.
         * @details Returns a string containing information about the matrix's dimensions, non-zeros, and properties.
         * @param[in] header Optional header string to prepend to the information.
         * @return A string containing matrix information.
         */
        std::string info(const std::string& header = "") const;

        /**
         * @brief Create a deep copy of the sparse matrix.
         * @details Creates a new sparse matrix with its own memory allocation
         *          and copies all elements from this matrix to the new matrix.
         * @return A new sparse matrix containing a copy of this matrix's data.
         */
        XxMatrix<T_Int,T_Scalar> copy() const;

        /**
         * @brief Create a reference copy (shallow copy) of the sparse matrix.
         * @details Creates a new sparse matrix object that references the same memory
         *          as this matrix. Changes to either matrix will affect both.
         * @return A sparse matrix that shares memory with this matrix.
         */
        XxMatrix<T_Int,T_Scalar> rcopy();

        /**
         * @brief Create a guarded reference copy (shallow copy) of the sparse matrix.
         * @details Creates a guarded sparse matrix object that references the same memory
         *          as this matrix. The guard ensures the reference is read-only.
         * @return A guarded sparse matrix that shares memory with this matrix.
         */
        Guard<XxMatrix<T_Int,T_Scalar>> rcopy() const;

        /**
         * @brief Move the sparse matrix's resources.
         * @details Transfers ownership of the memory to a new sparse matrix object,
         *          leaving this matrix in an empty state.
         * @return A new sparse matrix containing this matrix's memory.
         */
        XxMatrix<T_Int,T_Scalar> move();

        /**
         * @brief Output the sparse matrix to a stream.
         * @details Writes the sparse matrix elements to the specified output stream with specified precision.
         * @param[in,out] os The output stream to write to.
         * @param[in] prec The output precision (default: 0 for default precision).
         */
        void toStream(std::ostream& os, std::streamsize prec = 0) const;

        /**
         * @brief Scale the sparse matrix in-place.
         * @details Multiplies all non-zero elements of the sparse matrix by a scalar value.
         * @param[in] val The scalar value to multiply by.
         */
        void iscale(T_Scalar val);

        /**
         * @brief Transpose the sparse matrix.
         * @details Returns a transposed view of the sparse matrix.
         * @return A virtual transpose expression.
         */
        alias::VirtualTrans_csr<T_Int,T_Scalar> transpose() const;

        /**
         * @brief Conjugate transpose the sparse matrix.
         * @details Returns a conjugate transposed view of the sparse matrix.
         * @return A virtual conjugate transpose expression.
         */
        alias::VirtualTrans_csr<T_Int,T_Scalar> ctranspose() const;

        /**
         * @brief Compute the complex conjugate.
         * @details Returns a virtual expression containing the complex conjugate of each non-zero element.
         * @return A virtual expression with conjugated elements.
         */
        alias::VirtualConj_csr<T_Int,T_Scalar> conjugate() const;

        /**
         * @brief Conjugate the sparse matrix in-place.
         * @details Replaces all non-zero elements with their complex conjugates.
         */
        void iconjugate();

        /**
         * @brief Compute the 1-norm.
         * @details Computes the maximum absolute column sum.
         * @return The 1-norm of the sparse matrix.
         */
        T_RScalar normOne() const;

        /**
         * @brief Compute the infinity norm.
         * @details Computes the maximum absolute row sum.
         * @return The infinity norm of the sparse matrix.
         */
        T_RScalar normInf() const;

        /**
         * @brief Compute the maximum norm.
         * @details Computes the maximum absolute value of all non-zero elements.
         * @return The maximum norm of the sparse matrix.
         */
        T_RScalar normMax() const;

        /**
         * @brief Compute the Frobenius norm.
         * @details Computes the square root of the sum of squared absolute values of all non-zero elements.
         * @return The Frobenius norm of the sparse matrix.
         */
        T_RScalar normFro() const;

        /**
         * @brief Convert to general sparse matrix.
         * @details Creates a new general sparse matrix by copying all relevant elements from a structured sparse matrix.
         * @return A general sparse matrix containing the full data.
         */
        XxMatrix<T_Int,T_Scalar> general() const;

        /**
         * @brief Convert to dense matrix format.
         * @details Converts the sparse matrix to a dense matrix by expanding all elements.
         * @return A dense matrix containing the same data.
         */
        dns::XxMatrix<T_Scalar> toDns() const;

        /**
         * @brief Permute rows and columns.
         * @details Applies left and right permutations to the sparse matrix.
         * @param[in] P The row permutation matrix.
         * @param[in] Q The column permutation matrix.
         * @return A new sparse matrix containing the permuted elements.
         */
        XxMatrix<T_Int,T_Scalar> permuteLeftRight(const prm::PxMatrix<T_Int>& P, 
                                                  const prm::PxMatrix<T_Int>& Q) const;

        /**
         * @brief Permute rows.
         * @details Applies a left (row) permutation to the sparse matrix.
         * @param[in] P The row permutation matrix.
         * @return A new sparse matrix containing the permuted elements.
         */
        XxMatrix<T_Int,T_Scalar> permuteLeft(const prm::PxMatrix<T_Int>& P) const;

        /**
         * @brief Permute columns.
         * @details Applies a right (column) permutation to the sparse matrix.
         * @param[in] Q The column permutation matrix.
         * @return A new sparse matrix containing the permuted elements.
         */
        XxMatrix<T_Int,T_Scalar> permuteRight(const prm::PxMatrix<T_Int>& Q) const;

        /**
         * @brief Permute rows and columns symmetrically.
         * @details Applies the same permutation to both rows and columns (mirror permutation).
         * @param[in] P The permutation matrix to apply to both rows and columns.
         * @return A new sparse matrix containing the permuted elements.
         */
        XxMatrix<T_Int,T_Scalar> permuteMirror(const prm::PxMatrix<T_Int>& P) const;

        /**
         * @brief Extract a block as a new sparse matrix.
         * @details Creates a new sparse matrix containing a copy of a contiguous block of elements.
         * @param[in] ibgn The starting row index of the block.
         * @param[in] jbgn The starting column index of the block.
         * @param[in] ni The number of rows in the block.
         * @param[in] nj The number of columns in the block.
         * @return A new sparse matrix containing the block's data.
         */
        XxMatrix<T_Int,T_Scalar> block(T_Int ibgn, T_Int jbgn, T_Int ni, T_Int nj) const;

        /** @} */

        /** 
         * @name Creators/Generators
         * @{
         */

        /**
         * @brief Create a random sparse matrix.
         * @details Creates a sparse matrix with random non-zero values uniformly distributed in the specified range.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] nz The number of non-zero elements.
         * @param[in] pr The matrix property (default: General).
         * @param[in] lo The lower bound of the random values (default: 0).
         * @param[in] hi The upper bound of the random values (default: 1).
         * @return A sparse matrix containing random values.
         */
        static XxMatrix<T_Int,T_Scalar> random(T_Int nr, T_Int nc, T_Int nz, 
                                               const Property& pr = Property::General(),
                                               T_RScalar lo = T_RScalar(0), T_RScalar hi = T_RScalar(1));

        /**
         * @brief Create a view of existing CSR arrays.
         * @details Creates a guarded sparse matrix that references existing CSR format memory
         *          without taking ownership. The memory must remain valid for the lifetime
         *          of the returned view.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] rptr Pointer to the row pointer array.
         * @param[in] cidx Pointer to the column index array.
         * @param[in] vals Pointer to the values array.
         * @param[in] pr The matrix property (default: General).
         * @return A guarded sparse matrix that views the specified CSR arrays.
         */
        static Guard<XxMatrix<T_Int,T_Scalar>> view(T_Int nr, T_Int nc, 
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

        XxMatrix<T_Int,T_Scalar>& evaluateFrom(const alias::VirtualStrided_csr<T_Int,T_Scalar>& v)
        {
            if(*this) {
                v.evaluateOnExisting(*this);
            } else {
                v.evaluateOnNew(*this);
            }
            return *this;
        }

        XxMatrix<T_Int,T_Scalar>& evaluateFrom(const alias::GuardedStrided_csr<T_Int,T_Scalar>& gv)
        {
            return evaluateFrom(gv.get());
        }
};

/*-------------------------------------------------*/
} // namespace csr
} // namespace cla3p
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_stream_operators
 * @brief Writes to os the contents of mat.
 */
template <typename T_Int, typename T_Scalar>
std::ostream& operator<<(std::ostream& os, const cla3p::csr::XxMatrix<T_Int,T_Scalar>& mat)
{
    mat.toStream(os);
    return os;
}

/*-------------------------------------------------*/

#endif // CLA3P_CSR_XXMATRIX_HPP_