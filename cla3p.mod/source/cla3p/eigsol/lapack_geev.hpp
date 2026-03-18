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

#ifndef CLA3P_LAPACK_GEEV_HPP_
#define CLA3P_LAPACK_GEEV_HPP_

/**
 * @file
 */

#include <cla3p/types/integer.hpp>
#include <cla3p/support/heap_buffer.hpp>
#include <cla3p/dense/dns_cxvector.hpp>
#include <cla3p/dense/dns_cxmatrix.hpp>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Eigenvalue decomposition solver using LAPACK GEEV routines.
 * @details This class provides a high-level interface for computing eigenvalues and eigenvectors
 *          of general (non-symmetric) square matrices using the LAPACK GEEV implementation
 *          (real: `dgeev`/`sgeev`, complex: `zgeev`/`cgeev`).
 *
 *          The GEEV routine computes the eigenvalue decomposition of a general matrix @f$ A @f$:
 *          - Right eigenvectors: @f$ A v_i = \lambda_i v_i @f$
 *          - Left eigenvectors:  @f$ w_i^H A = \lambda_i w_i^H @f$
 *
 *          Results are always returned as complex matrices/vectors, even for real input matrices,
 *          since real matrices may have complex conjugate eigenvalue pairs. Internal buffers are
 *          managed automatically and can be pre-allocated via @ref reserve() to avoid repeated
 *          heap allocations when the solver is reused for multiple decompositions.
 *
 * @tparam T_Matrix The dense matrix type to decompose (can be real or complex).
 *
 * @section geev_usage Usage Example
 * @code
 * cla3p::dns::RdMatrix A = ...; // Input matrix
 *
 * LapackGeev<cla3p::dns::RdMatrix> geev;
 * geev.reserve(A.nrows()); // Optional: pre-allocate buffers
 * geev.decompose(A);
 *
 * const cla3p::dns::CdVector& eigenvalues      = geev.eigenvalues();
 * const cla3p::dns::CdMatrix& rightEigenvectors = geev.rightEigenvectors();
 * @endcode
 */
template <typename T_Matrix>
class LapackGeev {

    private:
        using T_Scalar = typename T_Matrix::value_type;
        using T_CScalar = typename TypeTraits<T_Scalar>::complex_type;
        using T_CVector = dns::CxVector<T_CScalar>;
        using T_CMatrix = dns::CxMatrix<T_CScalar>;

    public:

        // do not copy
        LapackGeev(const LapackGeev&) = delete;
        LapackGeev& operator=(const LapackGeev&) = delete;

        /**
         * @brief Default constructor.
         * @details Initializes the eigenvalue solver with default settings.
         *          @p calcLeft is set to @c false and @p calcRight is set to @c true,
         *          so only right eigenvectors are computed by default.
         */
        LapackGeev();

        /**
         * @brief Constructor with explicit eigenvector options.
         * @details Initializes the eigenvalue solver with the specified eigenvector
         *          computation options. Use this form when left eigenvectors are needed
         *          or when neither set of eigenvectors is required (eigenvalues only).
         * @param[in] calcLeft  If @c true, left eigenvectors will be computed in @ref decompose().
         * @param[in] calcRight If @c true, right eigenvectors will be computed in @ref decompose().
         */
        LapackGeev(bool calcLeft, bool calcRight);

        /**
         * @brief Destructor.
         * @details Releases all allocated resources and clears internal state.
         */
        ~LapackGeev();

        /**
         * @brief Pre-allocates internal buffers for a matrix of dimension @p n.
         * @details Allocates all scalar and complex workspace buffers required for
         *          eigenvalue decomposition of an @p n x @p n matrix. Calling this
         *          before @ref decompose() avoids repeated heap allocations when the
         *          same solver instance is reused across multiple decompositions of
         *          equal or smaller size.
         * @param[in] n Matrix dimension.
         */
        void reserve(int_t n);

        /**
         * @brief Clears all internal data and releases allocated memory.
         * @details Resets the solver to its initial state, releasing all computed results,
         *          internal workspace buffers, and resetting options to their defaults
         *          (@p calcLeft = @c false, @p calcRight = @c true).
         */
        void clear();

        /**
         * @brief Performs the eigenvalue decomposition of @p mat.
         * @details Computes the eigenvalues and, depending on the current settings,
         *          the left and/or right eigenvectors of @p mat using the LAPACK GEEV
         *          routine. The input matrix is copied internally before factorization;
         *          the original is left unmodified.
         *
         *          Results are accessible via @ref eigenvalues(), @ref leftEigenvectors(),
         *          and @ref rightEigenvectors() after this call returns.
         *
         * @param[in] mat Square general dense matrix to decompose.
         * @note For real input matrices, eigenvalues and eigenvectors are always returned as complex.
         */
        void decompose(const T_Matrix& mat);

        /**
         * @brief Sets whether to compute left eigenvectors.
         * @details Takes effect on the next call to @ref decompose().
         *          Changing this setting after @ref reserve() has been called will require
         *          new buffer allocations on the next @ref decompose().
         * @param[in] calcLeft If @c true, left eigenvectors will be computed.
         */
        void setCalcLeft(bool calcLeft);

        /**
         * @brief Sets whether to compute right eigenvectors.
         * @details Takes effect on the next call to @ref decompose().
         *          Changing this setting after @ref reserve() has been called will require
         *          new buffer allocations on the next @ref decompose().
         * @param[in] calcRight If @c true, right eigenvectors will be computed.
         */
        void setCalcRight(bool calcRight);

        /**
         * @brief Checks if left eigenvectors will be computed.
         * @return @c true if left eigenvectors will be computed, @c false otherwise.
         */
        bool getCalcLeft() const;

        /**
         * @brief Checks if right eigenvectors will be computed.
         * @return @c true if right eigenvectors will be computed, @c false otherwise.
         */
        bool getCalcRight() const;

        /**
         * @brief Returns the computed eigenvalues.
         * @details Each element @p i of the returned vector corresponds to the @p i-th
         *          eigenvalue @f$ \lambda_i @f$, matching the column ordering of
         *          @ref leftEigenvectors() and @ref rightEigenvectors().
         * @return Const reference to a complex vector of length @p n containing the eigenvalues.
         * @note Only valid after a successful call to @ref decompose().
         */
        const T_CVector& eigenvalues() const;

        /**
         * @brief Returns the computed left eigenvectors.
         * @details Column @p i of the returned matrix contains the left eigenvector
         *          @f$ w_i @f$ satisfying @f$ w_i^H A = \lambda_i w_i^H @f$, where
         *          @f$ \lambda_i @f$ is @ref eigenvalues() @p [i].
         * @return Const reference to a complex @p n x @p n matrix of left eigenvectors.
         * @note Only valid after calling @ref decompose() with @ref getCalcLeft() == @c true.
         */
        const T_CMatrix& leftEigenvectors() const;

        /**
         * @brief Returns the computed right eigenvectors.
         * @details Column @p i of the returned matrix contains the right eigenvector
         *          @f$ v_i @f$ satisfying @f$ A v_i = \lambda_i v_i @f$, where
         *          @f$ \lambda_i @f$ is @ref eigenvalues() @p [i].
         * @return Const reference to a complex @p n x @p n matrix of right eigenvectors.
         * @note Only valid after calling @ref decompose() with @ref getCalcRight() == @c true.
         */
        const T_CMatrix& rightEigenvectors() const;

    private:
        bool m_calcLeft;
        bool m_calcRight;

        T_Matrix m_matrixBackup;

        T_CVector m_eigenvalues;
        T_CMatrix m_leftEigenvectors;
        T_CMatrix m_rightEigenvectors;

        HeapBuffer<T_Scalar> m_scalarBuffers;
        HeapBuffer<T_CScalar> m_complexBuffers;

        void defaults();
        void clearInternalWrappers();
        void clearInternalBuffers();
        void resizeInternalObjects(int_t n);
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_LAPACK_GEEV_HPP_
