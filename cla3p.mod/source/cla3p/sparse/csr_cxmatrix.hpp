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

#ifndef CLA3P_CSR_CXMATRIX_HPP_
#define CLA3P_CSR_CXMATRIX_HPP_

/**
 * @file
 */

#include <cla3p/sparse/csr_xxmatrix.hpp>
#include <cla3p/virtuals/virtual_strided.hpp>

/*-------------------------------------------------*/
namespace cla3p { 
namespace csr {
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief The sparse complex matrix class (compressed sparse row format).
 * @details Represents a sparse complex-valued matrix in compressed sparse row (CSR) format.
 * @tparam T_Int The integer type for indices.
 * @tparam T_Scalar The complex scalar type (e.g., std::complex<float>, std::complex<double>).
 */
template <typename T_Int, typename T_Scalar>
class CxMatrix : public XxMatrix<T_Int,T_Scalar> {

    private:
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

    public:

        //
        // Convertors
        // Move convertors intentionally left as non-explicit
        //
        explicit CxMatrix(const XxMatrix<T_Int,T_Scalar>& other);
        CxMatrix<T_Int,T_Scalar>& operator=(const XxMatrix<T_Int,T_Scalar>& other);

        CxMatrix(XxMatrix<T_Int,T_Scalar>&& other);
        CxMatrix<T_Int,T_Scalar>& operator=(XxMatrix<T_Int,T_Scalar>&& other);

        template <typename T_Virtual>
        CxMatrix(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v) : XxMatrix<T_Int,T_Scalar>(v) {}
        template <typename T_Virtual>
        CxMatrix<T_Int,T_Scalar>& operator=(const alias::VirtualExpr_csr<T_Int,T_Scalar,T_Virtual>& v) { XxMatrix<T_Int,T_Scalar>::operator=(v); return *this; }

        /**
         * @name Constructors
         * @{
         */

        /**
         * @brief Default constructor.
         * @details Creates an empty complex sparse matrix with no allocated memory.
         */
        CxMatrix();

        /**
         * @brief Dimension constructor.
         * @details Creates a complex sparse matrix of the specified dimensions in CSR format and allocates memory.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] nz The number of non-zero elements.
         * @param[in] pr The matrix property (default: General).
         */
        explicit CxMatrix(T_Int nr, T_Int nc, T_Int nz, const Property& pr = Property::General());

        /**
         * @brief Auxiliary constructor.
         * @details Creates a complex sparse matrix using existing CSR format arrays.
         * @param[in] nr The number of rows.
         * @param[in] nc The number of columns.
         * @param[in] rptr Pointer to existing row pointer array.
         * @param[in] cidx Pointer to existing column index array.
         * @param[in] vals Pointer to existing values array.
         * @param[in] bind If true, the matrix takes ownership of the memory.
         * @param[in] pr The matrix property (default: General).
         */
        explicit CxMatrix(T_Int nr, T_Int nc, T_Int *rptr, T_Int *cidx, T_Scalar *vals, bool bind, const Property& pr = Property::General());

        /**
         * @brief Copy constructor.
         * @details Creates a new complex sparse matrix by copying another complex sparse matrix.
         * @param[in] other The complex sparse matrix to copy.
         */
        CxMatrix(const CxMatrix<T_Int,T_Scalar>& other) = default;

        /**
         * @brief Move constructor.
         * @details Creates a new complex sparse matrix by moving resources from another complex sparse matrix.
         * @param[in] other The complex sparse matrix to move from.
         */
        CxMatrix(CxMatrix<T_Int,T_Scalar>&& other) = default;

        /**
         * @brief Destructor.
         * @details Destroys the complex sparse matrix and releases allocated memory.
         */
        ~CxMatrix();

        /** @} */

        /**
         * @name Operators
         * @{
         */

        /**
         * @brief Copy assignment operator.
         * @details Copies the contents of another complex sparse matrix to this complex sparse matrix.
         * @param[in] other The complex sparse matrix to copy.
         * @return Reference to this complex sparse matrix.
         */
        CxMatrix<T_Int,T_Scalar>& operator=(const CxMatrix<T_Int,T_Scalar>& other) = default;

        /**
         * @brief Move assignment operator.
         * @details Moves resources from another complex sparse matrix to this complex sparse matrix.
         * @param[in] other The complex sparse matrix to move from.
         * @return Reference to this complex sparse matrix.
         */
        CxMatrix<T_Int,T_Scalar>& operator=(CxMatrix<T_Int,T_Scalar>&& other) = default;

        /** @} */

        /** 
         * @name Public Member Functions
         * @{
         */

        /**
         * @brief Accesses the real part of the complex sparse matrix.
         * @details Returns a strided view of the real components that can be used for both reading and writing.
         * @return A virtual strided view of the real part.
         */
        alias::VirtualStrided_csr<T_Int,T_RScalar> real();

        /**
         * @brief Accesses the imaginary part of the complex sparse matrix.
         * @details Returns a strided view of the imaginary components that can be used for both reading and writing.
         * @return A virtual strided view of the imaginary part.
         */
        alias::VirtualStrided_csr<T_Int,T_RScalar> imag();

        /**
         * @brief Accesses the real part of the complex sparse matrix (const).
         * @details Returns a read-only strided view of the real components.
         * @return A guarded strided view of the real part.
         */
        alias::GuardedStrided_csr<T_Int,T_RScalar> real() const;

        /**
         * @brief Accesses the imaginary part of the complex sparse matrix (const).
         * @details Returns a read-only strided view of the imaginary components.
         * @return A guarded strided view of the imaginary part.
         */
        alias::GuardedStrided_csr<T_Int,T_RScalar> imag() const;

        /** @} */
};

/*-------------------------------------------------*/
} // namespace csr
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSR_CXMATRIX_HPP_
