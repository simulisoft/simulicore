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

#ifndef CLA3P_LAPACK_SYEV_HPP_
#define CLA3P_LAPACK_SYEV_HPP_

/**
 * @file
 */

#include <cla3p/types/integer.hpp>
#include <cla3p/support/heap_buffer.hpp>
#include <cla3p/dense/dns_xxvector.hpp>

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Eigenvalue decomposition solver using LAPACK SYEV/HEEV routines.
 * @details This class provides a high-level interface for computing eigenvalues and
 *          eigenvectors of symmetric (real) or Hermitian (complex) square matrices
 *          using the LAPACK HEEV implementation (real: `dsyev`/`ssyev`, complex: `zheev`/`cheev`).
 *
 *          The HEEV routine computes the eigenvalue decomposition of a symmetric/Hermitian
 *          matrix @f$ A @f$:
 *          @f[
 *            A \, v_i = \lambda_i \, v_i, \qquad \lambda_i \in \mathbb{R}
 *          @f]
 *          where all eigenvalues @f$ \lambda_i @f$ are guaranteed to be real. The
 *          eigenvectors form an orthonormal basis satisfying @f$ V^H V = I @f$.
 *
 *          Unlike @ref LapackGeev, eigenvalues are returned as real-valued quantities.
 *          Eigenvectors share the scalar type of the input matrix (real for real input,
 *          complex for complex input). Internal buffers are managed automatically and can
 *          be pre-allocated via @ref reserve() to avoid repeated heap allocations when
 *          the solver is reused for multiple decompositions.
 *
 * @tparam T_Matrix The dense symmetric/Hermitian matrix type to decompose.
 *
 * @section syev_usage Usage Example
 * @code
 * cla3p::dns::RdMatrix A = ...; // Input symmetric matrix
 *
 * LapackSyev<cla3p::dns::RdMatrix> syev;
 * syev.reserve(A.nrows()); // Optional: pre-allocate buffers
 * syev.decompose(A);
 *
 * const cla3p::dns::RdVector& eigenvalues  = syev.eigenvalues();
 * const cla3p::dns::RdMatrix& eigenvectors = syev.eigenvectors();
 * @endcode
 */
template <typename T_Matrix>
class LapackSyev {

    private:
        using T_Scalar = typename T_Matrix::value_type;
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
        using T_RVector = dns::XxVector<T_RScalar>;

    public:

        // do not copy
        LapackSyev(const LapackSyev&) = delete;
        LapackSyev& operator=(const LapackSyev&) = delete;

        /**
         * @brief Default constructor.
         * @details Initializes the eigenvalue solver with default settings.
         *          @p calcEigenvectors is set to @c true, so eigenvectors are computed by default.
         */
        LapackSyev();

        /**
         * @brief Constructor with explicit eigenvector option.
         * @details Initializes the eigenvalue solver with the specified eigenvector
         *          computation option. Use this form when only eigenvalues are needed,
         *          which avoids the additional cost of computing the eigenvector matrix.
         * @param[in] calcEigenvectors If @c true, eigenvectors will be computed in @ref decompose().
         */
        LapackSyev(bool calcEigenvectors);

        /**
         * @brief Destructor.
         * @details Releases all allocated resources and clears internal state.
         */
        ~LapackSyev();

        /**
         * @brief Clears all internal data and releases allocated memory.
         * @details Resets the solver to its initial state, releasing all computed results,
         *          internal workspace buffers, and resetting options to their defaults
         *          (@p calcEigenvectors = @c true).
         */
        void clear();

        /**
         * @brief Pre-allocates internal buffers for a matrix of dimension @p n.
         * @details Allocates all real and scalar workspace buffers required for
         *          eigenvalue decomposition of an @p n x @p n matrix. Calling this
         *          before @ref decompose() avoids repeated heap allocations when the
         *          same solver instance is reused across multiple decompositions of
         *          equal or smaller size.
         * @param[in] n Matrix dimension.
         */
        void reserve(int_t n);

        /**
         * @brief Performs the eigenvalue decomposition of @p mat.
         * @details Computes the eigenvalues and, depending on the current setting,
         *          the eigenvectors of @p mat using the LAPACK HEEV routine. The input
         *          matrix is copied internally before factorization; the original is
         *          left unmodified.
         *
         *          Results are accessible via @ref eigenvalues() and @ref eigenvectors()
         *          after this call returns.
         *
         * @param[in] mat Square symmetric/Hermitian dense matrix to decompose.
         * @note Eigenvalues are always real, regardless of whether the input is real or complex.
         */
        void decompose(const T_Matrix& mat);

        /**
         * @brief Sets whether to compute eigenvectors.
         * @details Takes effect on the next call to @ref decompose().
         *          Changing this setting after @ref reserve() has been called will require
         *          new buffer allocations on the next @ref decompose().
         * @param[in] calcEigenvectors If @c true, eigenvectors will be computed.
         */
        void setCalcEigenvectors(bool calcEigenvectors);

        /**
         * @brief Checks if eigenvectors will be computed.
         * @return @c true if eigenvectors will be computed, @c false otherwise.
         */
        bool getCalcEigenvectors() const;

        /**
         * @brief Returns the computed eigenvalues.
         * @details Each element @p i of the returned vector is the @p i-th eigenvalue
         *          @f$ \lambda_i @f$, sorted in ascending order, matching the column
         *          ordering of @ref eigenvectors().
         * @return Const reference to a real vector of length @p n containing the eigenvalues.
         * @note Only valid after a successful call to @ref decompose().
         */
        const T_RVector& eigenvalues() const;

        /**
         * @brief Returns the computed eigenvectors.
         * @details Column @p i of the returned matrix contains the orthonormal eigenvector
         *          @f$ v_i @f$ satisfying @f$ A v_i = \lambda_i v_i @f$, where
         *          @f$ \lambda_i @f$ is @ref eigenvalues() @p [i]. The columns form an
         *          orthonormal basis: @f$ V^H V = I @f$.
         * @return Const reference to an @p n x @p n matrix of eigenvectors, with the same
         *         scalar type as the input matrix.
         * @note Only valid after calling @ref decompose() with @ref getCalcEigenvectors() == @c true.
         */
        const T_Matrix& eigenvectors() const;

    private:
        bool m_calcEigenvectors;

        T_Matrix  m_matrixBackup;
        T_RVector m_eigenvalues;
        T_Matrix  m_eigenvectors;

        HeapBuffer<T_RScalar> m_realBuffers;
        HeapBuffer<T_Scalar> m_scalarBuffers;

        void clearInternalWrappers();
        void clearInternalBuffers();
        void resizeInternalObjects(int_t n);

        void defaults();
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_LAPACK_SYEV_HPP_
