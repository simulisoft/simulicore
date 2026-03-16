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

#ifndef CULITE_DEFAULT_SVD_HPP_
#define CULITE_DEFAULT_SVD_HPP_

/**
 * @file
 */

#include <culite/support/imalloc.hpp>
#include <culite/generic/cusolver_handler.hpp>
#include <culite/dense/dns_xxvector.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief Singular Value Decomposition (SVD) solver using cuSOLVER.
 * @details This class provides a high-level interface for computing the singular value decomposition
 *          of general matrices using the cuSOLVER library.
 *          
 *          For a given @f$ m \times n @f$ matrix A, the SVD computes the factorization:
 *          @f[
 *              A = U \Sigma V^T
 *          @f]
 *          where:
 *          - @f$ U @f$ is an @f$ m \times m @f$ (or @f$ m \times \min(m,n) @f$ if limited) orthogonal matrix of left singular vectors
 *          - @f$ \Sigma @f$ is an @f$ m \times n @f$ diagonal matrix with non-negative real singular values
 *          - @f$ V^T @f$ is the transpose of an @f$ n \times n @f$ (or @f$ \min(m,n) \times n @f$ if limited) orthogonal matrix
 * 
 *          The computation of singular vectors can be controlled independently via policy settings
 *          (Full, Limited, or None) to optimize performance and memory usage.
 * 
 * @tparam T_Matrix The matrix type (must be a dense GPU matrix type).
 * 
 * @section svd_usage Usage Example
 * @code
 * culite::dns::RdMatrix A = ...; // Input matrix
 * 
 * DefaultSVD<culite::dns::RdMatrix> svd;
 * svd.reserve(A); // Optional
 * svd.decompose(A);
 * 
 * const culite::dns::RdVector& S = svd.singularValues();
 * const culite::dns::RdMatrix& U = svd.leftSingularVectors();
 * const culite::dns::RdMatrix& V = svd.rightSingularVectors();
 * @endcode
 */
template <typename T_Matrix>
class DefaultSVD {

    private:
        using T_Scalar = typename T_Matrix::value_type;
        using T_RScalar = typename TypeTraits<T_Scalar>::real_type;
        using T_RVector = dns::XxVector<T_RScalar>;

    public:

        // do not copy
        DefaultSVD(const DefaultSVD&) = delete;
        DefaultSVD& operator=(const DefaultSVD&) = delete;

        /**
         * @brief Default constructor.
         * @details Initializes the SVD solver with default settings.
         *          Both leftPolicy and rightPolicy are initialized to Limited.
         * @param[in] cusolver Reference to a cuSOLVER handler instance (defaults to global handler).
         */
        DefaultSVD(CuSolverHandler& cusolver = globalCuSolverHandler());

        /**
         * @brief Constructor with specified policies.
         * @details Initializes the SVD solver with the specified cuSOLVER handler and
         *          singular vector computation policies.
         * @param[in] leftPolicy Policy for computing left singular vectors U (Limited, Full, or None).
         * @param[in] rightPolicy Policy for computing right singular vectors V (Limited, Full, or None).
         * @param[in] cusolver Reference to a cuSOLVER handler instance (defaults to global handler).
         */
        DefaultSVD(svdPolicy_t leftPolicy, 
                   svdPolicy_t rightPolicy, 
                   CuSolverHandler& cusolver = globalCuSolverHandler());
        
        /**
         * @brief Destructor.
         * @details Releases all allocated resources and clears internal state.
         */
        ~DefaultSVD();

        /**
         * @brief Reserves workspace memory for SVD decomposition.
         * @details Allocates internal buffers required for SVD computation based on
         *          the dimensions of the input matrix and current policies. Must be called before @ref decompose.
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
         * @brief Gets the current policy for computing left singular vectors.
         * @return The current left singular vector computation policy.
         */
        svdPolicy_t getLeftPolicy () const { return m_LeftPolicy; }
        
        /**
         * @brief Gets the current policy for computing right singular vectors.
         * @return The current right singular vector computation policy.
         */
        svdPolicy_t getRightPolicy() const { return m_RightPolicy; }

        /**
         * @brief Sets the policy for computing left singular vectors.
         * @details Changing the policy requires calling @ref reserve again before the next decomposition.
         * @param[in] policy The desired left singular vector computation policy.
         */
        void setLeftPolicy (svdPolicy_t policy) { m_LeftPolicy = policy; }
        
        /**
         * @brief Sets the policy for computing right singular vectors.
         * @details Changing the policy requires calling @ref reserve again before the next decomposition.
         * @param[in] policy The desired right singular vector computation policy.
         */
        void setRightPolicy(svdPolicy_t policy) { m_RightPolicy = policy; }

        /**
         * @brief Performs singular value decomposition on the input matrix.
         * @details Computes the singular values and optionally the left and/or right singular vectors
         *          of the input matrix using cuSOLVER's SVD routine.
         *          
         *          The input matrix is not modified. Results can be retrieved using the
         *          accessor methods: @ref singularValues, @ref leftSingularVectors, @ref rightSingularVectors.
         * 
         * @param[in] mat The matrix to decompose.
         * 
         * @pre @ref reserve must have been called with a matrix of matching dimensions.
         * @note Singular values are returned in descending order.
         */
        void decompose(const T_Matrix& mat);

        /**
         * @brief Gets the computed singular values.
         * @return Const reference to the vector containing the singular values in descending order.
         * @note Only valid after calling @ref decompose.
         * @note Singular values satisfy @f$ \sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_k \geq 0 @f$ where @f$ k = \min(m,n) @f$.
         */
        const T_RVector& singularValues() const { return m_singularValues; }
        
        /**
         * @brief Gets the computed left singular vectors.
         * @return Const reference to the matrix containing the left singular vectors (U matrix).
         * @note Only valid after calling @ref decompose.
         * @note Each column i contains the left singular vector corresponding to singular value i.
         * @note Dimensions depend on left policy: Limited gives @f$ m \times \min(m,n) @f$, Full gives @f$ m \times m @f$.
         */
        const T_Matrix& leftSingularVectors() const { return m_leftSingularVectors; }
        
        /**
         * @brief Gets the computed right singular vectors.
         * @return Const reference to the matrix containing the right singular vectors (V matrix).
         * @note Only valid after calling @ref decompose.
         * @note Each column i contains the right singular vector corresponding to singular value i.
         * @note Returns V (not V^T). Dimensions depend on right policy: Limited gives @f$ n \times \min(m,n) @f$, Full gives @f$ n \times n @f$.
         */
        const T_Matrix& rightSingularVectors() const { return m_rightSingularVectors; }

    private:
        CuSolverHandler& m_cusolver;
        svdPolicy_t m_LeftPolicy;
        svdPolicy_t m_RightPolicy;

        T_RVector m_singularValues;
        T_Matrix m_leftSingularVectors;
        T_Matrix m_rightSingularVectors;

        DeviceBuffer<T_RScalar> m_deviceValuesBuffer;
        DeviceBuffer<T_Scalar> m_deviceVectorBuffer;

        void defaults();
        void clearOutput();
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_DEFAULT_SVD_HPP_
