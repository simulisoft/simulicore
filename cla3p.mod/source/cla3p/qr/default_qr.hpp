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

#ifndef CLA3P_DEFAULT_QR_HPP_
#define CLA3P_DEFAULT_QR_HPP_

/**
 * @file
 */

#include "cla3p/types/integer.hpp"
#include "cla3p/types/enums.hpp"
#include "cla3p/support/heap_buffer.hpp"
#include "cla3p/dense/dns_xxvector.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief QR Decomposition implementation for dense matrices.
 * @tparam T_Matrix The dense matrix type to decompose.
 * @details This class provides a complete QR decomposition interface,
 *          decomposing a matrix @f$ A @f$ into @f$ A = Q \cdot R @f$,
 *          where @f$ Q @f$ is an orthogonal (or unitary) matrix and @f$ R @f$ is upper triangular.
 */
template <typename T_Matrix>
class DefaultQR {

    private:
        using T_Scalar = typename T_Matrix::value_type;
        using T_Vector = dns::XxVector<T_Scalar>;

    public:

        /**
         * @brief Constructs an uninitialized QR decomposition object.
         * @details Initializes an empty QR object with no allocated memory for decomposition results.
         */
        DefaultQR();

        /**
         * @brief Constructs a QR object with pre-allocated buffers.
         * @details Initializes the QR decomposition object and allocates internal buffers
         *          to accommodate matrices up to the specified dimensions.
         * @param[in] m The maximum number of matrix rows.
         * @param[in] n The maximum number of matrix columns.
         * @param[in] qrPolicy The computation policy for the orthogonal matrix @f$ Q @f$.
         */
        DefaultQR(int_t m, int_t n, qrPolicy_t qrPolicy = qrPolicy_t::Reflection);

        /**
         * @brief Destroys the QR decomposition object.
         * @details Releases all allocated memory and clears internal data structures.
         */
        ~DefaultQR();

        /**
         * @brief Clears all decomposition results and resets the object.
         * @details Deallocates decomposition results and returns the object to its default state.
         */
        void clear();

        /**
         * @brief Retrieves the orthogonal matrix computation policy.
         * @details Returns the policy that determines whether matrix @f$ Q @f$ will be computed explicitly
         *          or represented implicitly as a product of elementary Householder reflectors.
         * @return The current orthogonal matrix computation policy.
         */
        qrPolicy_t getQrPolicy() const;

        /**
         * @brief Configures the orthogonal matrix computation policy.
         * @details Specifies whether matrix @f$ Q @f$ should be computed explicitly or
         *          represented implicitly as a product of elementary Householder reflectors.
         * @param[in] qrPolicy The orthogonal matrix computation policy.
         */
        void setQrPolicy(qrPolicy_t qrPolicy);

        /**
         * @brief Pre-allocates buffers for decomposition.
         * @details Allocates internal buffers to accommodate matrices with up to @p m rows and @p n columns.
         *          This avoids memory reallocation during subsequent decompositions.
         * @param[in] m The maximum number of matrix rows to support.
         * @param[in] n The maximum number of matrix columns to support.
         */
        void reserve(int_t m, int_t n);

        /**
         * @brief Performs QR Decomposition on the input matrix.
         * @details Computes the decomposition @f$ A = Q \cdot R @f$ of the input matrix.
         *          The computation respects the policy set via @ref setQrPolicy().
         * @param[in] mat The matrix to decompose.
         */
        void decompose(const T_Matrix& mat);

        /**
         * @brief Retrieves the upper triangular matrix from the decomposition.
         * @details Returns the matrix @f$ R @f$ of size @f$ \min(m,n) \times n @f$
         *          from the QR decomposition @f$ A = Q \cdot R @f$.
         * @return A constant reference to the upper triangular matrix @f$ R @f$.
         */
        const T_Matrix& R() const;

        /**
         * @brief Retrieves the orthogonal matrix from the decomposition.
         * @details Returns the matrix @f$ Q @f$ of size @f$ m \times \min(m,n) @f$
         *          from the QR decomposition @f$ A = Q \cdot R @f$.
         *          Only available when QR policy is set to Full; otherwise returns a reference to an empty matrix.
         * @return A constant reference to the orthogonal matrix @f$ Q @f$.
         */
        const T_Matrix& Q() const;

        /**
         * @brief Retrieves the Householder reflection coefficients.
         * @details Returns the vector @f$ \tau @f$ of size @f$ \min(m,n) @f$ containing
         *          the scalar factors for elementary Householder reflectors used to construct @f$ Q @f$.
         * @return A constant reference to the Householder scalar coefficient vector.
         */
        const T_Vector& tau() const;

        /**
         * @brief Retrieves the matrix of elementary Householder reflectors.
         * @details Returns a matrix containing the elementary Householder vectors used to implicitly represent @f$ Q @f$.
         *          This is available when the QR policy is Reflection.
         * @return A constant reference to the elementary reflector matrix.
         */
        const T_Matrix& elementaryReflectors() const;

    private:
        qrPolicy_t m_qrPolicy;

        T_Matrix m_matrixA;
        T_Vector m_tauVector;
        T_Matrix m_matrixR;
        T_Matrix m_matrixQ;

        HeapBuffer<T_Scalar> m_scalarBuffers;

        void defaults();

        void clearInternalWrappers();
        void clearInternalBuffers();
        void resizeInternalObjects(const T_Matrix& mat);

        void fillMatrixR();
        void fillMatrixQ();
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_DEFAULT_QR_HPP_
