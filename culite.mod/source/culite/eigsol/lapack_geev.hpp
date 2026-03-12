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

#ifndef CULITE_LAPACK_GEEV_HPP_
#define CULITE_LAPACK_GEEV_HPP_

/**
 * @file
 */

#include <culite/support/imalloc.hpp>
#include <culite/generic/cusolver_handler.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Eigenvalue decomposition solver using LAPACK/cuSOLVER GEEV routines.
 * @details This class provides a high-level interface for computing eigenvalues and eigenvectors
 *          of general (non-symmetric) matrices using the cuSOLVER library's GEEV implementation.
 *          
 *          The GEEV routine computes the eigenvalue decomposition of a general matrix A:
 *          - Right eigenvectors: @f$ A v_i = \lambda_i v_i @f$
 *          - Left eigenvectors: @f$ w_i^H A = \lambda_i w_i^H @f$
 * 
 *          Results are always returned as complex matrices/vectors, even for real input matrices.
 * 
 * @tparam T_Matrix The matrix type (can be real or complex).
 * 
 * @section geev_usage Usage Example
 * @code
 * CuSolverHandler cusolver;
 * culite::dns::RdMatrix A = ...; // Input matrix
 * 
 * LapackGeev<culite::dns::RdMatrix> geev(cusolver, false, true);
 * geev.reserve(A);
 * geev.decompose(A);
 * 
 * auto eigenvalues = geev.eigenvalues();
 * auto rightEigenvectors = geev.rightEigenvectors();
 * @endcode
 */
template <typename T_Matrix>
class LapackGeev {

    private:
        using T_Scalar = typename T_Matrix::value_type;
        using T_CScalar = typename TypeTraits<T_Scalar>::complex_type;

    public:

        // do not copy
        LapackGeev(const LapackGeev&) = delete;
        LapackGeev& operator=(const LapackGeev&) = delete;

        /**
         * @brief Default constructor.
         * @details Initializes the eigenvalue solver with default settings.
         *          calcLeft is initialized to false and calcRight is initialized to true.
         * @param[in] cusolver Reference to a cuSOLVER handler instance (defaults to global handler).
         */
        LapackGeev(CuSolverHandler& cusolver = globalCuSolverHandler());

        /**
         * @brief Constructor.
         * @details Initializes the eigenvalue solver with the specified cuSOLVER handler
         *          and eigenvector computation options.
         * @param[in] calcLeft If true, computes left eigenvectors.
         * @param[in] calcRight If true, computes right eigenvectors.
         * @param[in] cusolver Reference to a cuSOLVER handler instance (defaults to global handler).
         */
        LapackGeev(bool calcLeft, bool calcRight, CuSolverHandler& cusolver = globalCuSolverHandler());
        
        /**
         * @brief Destructor.
         * @details Releases all allocated resources and clears internal state.
         */
        ~LapackGeev();

        /**
         * @brief Checks if left eigenvectors will be computed.
         * @return True if left eigenvectors will be computed, false otherwise.
         */
        bool getCalcLeft() const { return m_calcLeft; }
        
        /**
         * @brief Checks if right eigenvectors will be computed.
         * @return True if right eigenvectors will be computed, false otherwise.
         */
        bool getCalcRight() const { return m_calcRight; }

        /**
         * @brief Sets whether to compute left eigenvectors.
         * @param[in] calcLeft If true, left eigenvectors will be computed in @ref decompose.
         */
        void setCalcLeft(bool calcLeft) { m_calcLeft = calcLeft; }
        
        /**
         * @brief Sets whether to compute right eigenvectors.
         * @param[in] calcRight If true, right eigenvectors will be computed in @ref decompose.
         */
        void setCalcRight(bool calcRight) { m_calcRight = calcRight; }

        /**
         * @brief Reserves workspace memory for eigenvalue decomposition.
         * @details Allocates internal buffers required for eigenvalue computation based on
         *          the dimensions of the input matrix. Must be called before @ref decompose.
         * @param[in] mat The matrix for which to reserve workspace (determines dimensions).
         */
        void reserve(const T_Matrix& mat);

        /**
         * @brief Clears all internal data and releases allocated memory.
         * @details Resets the solver to its initial state, releasing all computed results
         *          and internal workspace buffers.
         */
        void clear();

        /**
         * @brief Performs eigenvalue decomposition on the input matrix.
         * @details Computes the eigenvalues and optionally the left and/or right eigenvectors
         *          of the input matrix using cuSOLVER's GEEV routine.
         *          
         *          The input matrix is not modified. Results can be retrieved using the
         *          accessor methods: @ref eigenvalues, @ref leftEigenvectors, @ref rightEigenvectors.
         * 
         * @param[in] mat The matrix to decompose.
         * 
         * @pre @ref reserve must have been called with a matrix of matching dimensions.
         * @note For real input matrices, eigenvalues and eigenvectors are returned as complex.
         */
        void decompose(const T_Matrix& mat);

        /**
         * @brief Gets the computed eigenvalues.
         * @return Const reference to the complex vector containing the eigenvalues.
         * @note Only valid after calling @ref decompose.
         */
        const dns::CxVector<T_CScalar>& eigenvalues() const { return m_eigenvalues; }
        
        /**
         * @brief Gets the computed left eigenvectors.
         * @return Const reference to the complex matrix containing the left eigenvectors.
         * @note Only valid after calling @ref decompose with calcLeft=true.
         * @note Each column i contains the left eigenvector corresponding to eigenvalue i.
         */
        const dns::CxMatrix<T_CScalar>& leftEigenvectors() const { return m_leftEigenvectors; }
        
        /**
         * @brief Gets the computed right eigenvectors.
         * @return Const reference to the complex matrix containing the right eigenvectors.
         * @note Only valid after calling @ref decompose with calcRight=true.
         * @note Each column i contains the right eigenvector corresponding to eigenvalue i.
         */
        const dns::CxMatrix<T_CScalar>& rightEigenvectors() const { return m_rightEigenvectors; }

    private:
        CuSolverHandler& m_cusolver;
        bool m_calcLeft;
        bool m_calcRight;

        dns::CxVector<T_CScalar> m_eigenvalues;
        dns::CxMatrix<T_CScalar> m_leftEigenvectors;
        dns::CxMatrix<T_CScalar> m_rightEigenvectors;

        DeviceBuffer<T_CScalar> m_deviceBuffer;

        void defaults();
        void clearOutput();
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_LAPACK_GEEV_HPP_
