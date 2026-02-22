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

#ifndef CCULITE_CUSOLVER_HANDLER_HPP_
#define CCULITE_CUSOLVER_HANDLER_HPP_

/**
 * @file
 */

#include <cusolverDn.h>
#include <cla3p/checks/basic_checks.hpp>

#include "culite/types/traits.hpp"
#include "culite/error/cuda.hpp"
#include "culite/support/imalloc.hpp"
#include "culite/support/utils.hpp"
#include "culite/bulk/dns1D.hpp"
#include "culite/bulk/dns2D.hpp"
#include "culite/dense/dns_cxvector.hpp"
#include "culite/dense/dns_cxmatrix.hpp"

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The cuSolver handler class.
 * @details This class provides a wrapper around the cuSOLVER library for performing
 *          linear algebra operations on GPU devices. It manages the cuSOLVER handle
 *          and internal workspace buffers required for factorization and solve operations.
 */
class CuSolverHandler {

    private:
        using cuSolverInt = int64_t;
        using infoInt = int;

    public:

        // disable copying
        CuSolverHandler(const CuSolverHandler&) = delete;
        CuSolverHandler& operator=(const CuSolverHandler&) = delete;

        /**
         * @brief Constructor.
         * @details Initializes the cuSOLVER handle and internal state.
         */
        CuSolverHandler();

        
        /**
         * @brief Destructor.
         * @details Destroys the cuSOLVER handle and releases all allocated resources.
         */
        ~CuSolverHandler();
        
        /**
         * @brief Get the cuSolver handle.
         * @return The cuSolver handle.
         */
        cusolverDnHandle_t handle() { return m_handle; }

        /**
         * @brief Clears all internal workspace buffers.
         * @details Releases memory allocated for pivot indices, info arrays, and workspace buffers,
         *          resetting the handler to its initial state.
         */
        void clear();

        /**
         * @brief Reserves workspace memory for LU factorization.
         * @details Computes the required workspace size and allocates buffers for performing
         *          LU decomposition on matrix @p A. This includes memory for pivot indices,
         *          info arrays, factorization storage, and device/host workspaces.
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix for which to reserve workspace.
         */
        template <typename T_Matrix>
        void reserveLU(const T_Matrix& A)
        {
            using T_Scalar = typename T_Matrix::value_type;
            cusolverDnParams_t params = nullptr;
            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrf_bufferSize(handle(),
                                        params,
                                        A.nrows(), A.ncols(),
                                        TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                        TypeTraits<T_Scalar>::cuda_type(),
                                        &m_workspaceInBytesOnDevice,
                                        &m_workspaceInBytesOnHost);
            err::check_cusolver(cusolverStatus);

            ipiv().reserve(std::min(A.nrows(), A.ncols()));
            info().reserve(1);
            customWork().reserve(A.nrows() * A.ncols() * sizeof(T_Scalar));
            deviceWork().reserve(m_workspaceInBytesOnDevice);
            hostWork().reserve(m_workspaceInBytesOnHost);
        }

        /**
         * @brief Performs LU factorization with partial pivoting.
         * @details Computes the LU decomposition @f$ P A = L U @f$ of matrix @p A, where
         *          @p P is a permutation matrix, @p L is lower triangular, and @p U is upper triangular.
         *          The factorization is stored internally for subsequent solve operations.
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix to factorize.
         */
        template <typename T_Matrix>
        void decomposeLU(const T_Matrix& A)
        {
            using T_Scalar = typename T_Matrix::value_type;
            cuSolverInt n = A.nrows();
            memCopyD2D<T_Scalar>(A.nrows(), A.ncols(), 
                                 A.values(), A.ld(), 
                                 static_cast<T_Scalar*>(customWork().data()), n);
            cusolverDnParams_t params = nullptr;
            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrf(handle(),
                             params,
                             A.nrows(), A.ncols(),
                             TypeTraits<T_Scalar>::cuda_type(), customWork().data(), n,
                             ipiv().data(),
                             TypeTraits<T_Scalar>::cuda_type(),
                             deviceWork().data(), m_workspaceInBytesOnDevice,
                             hostWork().data(), m_workspaceInBytesOnHost,
                             info().data());
            err::check_cusolver(cusolverStatus);

            m_problemCudaType = TypeTraits<T_Scalar>::cuda_type();
            m_problemDim = n;
        }

        /**
         * @brief Solves a linear system using the computed LU factorization.
         * @details Solves the system @f$ A X = B @f$ using the previously computed LU factorization.
         *          The solution is stored in-place in matrix @p B.
         * @tparam T_Matrix The matrix type.
         * @param[in,out] B On input, the right-hand side matrix; on output, the solution matrix.
         */
        template <typename T_Matrix>
        void solveLU(T_Matrix& B)
        {
            cuSolverInt n = B.nrows();
            ::cla3p::similarity_dim_check(m_problemDim, n);
            using T_Scalar = typename T_Matrix::value_type;
            cusolverDnParams_t params = nullptr;
            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrs(handle(),
                             params,
                             CUBLAS_OP_N,
                             n, B.ncols(),
                             m_problemCudaType, customWork().data(), n,
                             ipiv().data(),
                             TypeTraits<T_Scalar>::cuda_type(), B.values(), B.ld(),
                             info().data());
            err::check_cusolver(cusolverStatus);
        }

        /**
         * @brief Reserves workspace memory for eigenvalue decomposition.
         * @details Computes the required workspace size and allocates buffers for performing
         *          eigenvalue decomposition on matrix @p A. This includes memory for eigenvalues,
         *          eigenvectors (if requested), and device/host workspaces.
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix for which to reserve workspace.
         * @param[in] calcLeft If true, reserves memory for left eigenvectors.
         * @param[in] calcRight If true, reserves memory for right eigenvectors.
         */
        template <typename T_Matrix>
        void reserveGeev(const T_Matrix& A, bool calcLeft, bool calcRight)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_CScalar = typename TypeTraits<T_Scalar>::complex_type;
            cusolverDnParams_t params = nullptr;

            cuSolverInt n = A.nrows();
            cuSolverInt sizeA = n * n * sizeof(T_Scalar);
            //cuSolverInt sizeW = (TypeTraits<T_Scalar>::is_real() ? 2 * n : n) * sizeof(T_Scalar);
            cuSolverInt sizeW = n * sizeof(T_CScalar);
            cuSolverInt sizeVL = calcLeft ? sizeA : 0;
            cuSolverInt sizeVR = calcRight ? sizeA : 0;

            customWork().reserve(sizeA + sizeW + sizeVL + sizeVR);

            const void *W = customWork().data();
            const void *VL = calcLeft ? W + sizeW : nullptr;
            const void *VR = calcRight ? W + sizeW + sizeVL : nullptr;

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgeev_bufferSize(handle(),
                                       params,
                                       calcLeft ? cusolverEigMode_t::CUSOLVER_EIG_MODE_VECTOR : cusolverEigMode_t::CUSOLVER_EIG_MODE_NOVECTOR,
                                       calcRight ? cusolverEigMode_t::CUSOLVER_EIG_MODE_VECTOR : cusolverEigMode_t::CUSOLVER_EIG_MODE_NOVECTOR,
                                       n,
                                       TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                       TypeTraits<T_CScalar>::cuda_type(), W,
                                       TypeTraits<T_Scalar>::cuda_type(), VL, n,
                                       TypeTraits<T_Scalar>::cuda_type(), VR, n,
                                       TypeTraits<T_Scalar>::cuda_type(),
                                       &m_workspaceInBytesOnDevice,
                                       &m_workspaceInBytesOnHost);

            err::check_cusolver(cusolverStatus);

            info().reserve(1);
            deviceWork().reserve(m_workspaceInBytesOnDevice);
            hostWork().reserve(m_workspaceInBytesOnHost);
        }

        /**
         * @brief Computes eigenvalues and eigenvectors of a general matrix.
         * @details Performs eigenvalue decomposition on matrix @p A, computing eigenvalues
         *          and optionally left and/or right eigenvectors. For a matrix @f$ A @f$,
         *          computes eigenvalues @f$ \lambda @f$ and eigenvectors such that
         *          @f$ A v = \lambda v @f$ (right) or @f$ w^H A = \lambda w^H @f$ (left).
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix for which to compute eigenvalues/eigenvectors.
         * @param[in] calcLeft If true, computes left eigenvectors.
         * @param[in] calcRight If true, computes right eigenvectors.
         */
        template <typename T_Matrix>
        void executeGeev(const T_Matrix& A, bool calcLeft, bool calcRight)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_CScalar = typename TypeTraits<T_Scalar>::complex_type;
            cusolverDnParams_t params = nullptr;

            cuSolverInt n = A.nrows();
            cuSolverInt sizeA = n * n * sizeof(T_Scalar);
            //cuSolverInt sizeW = (TypeTraits<T_Scalar>::is_real() ? 2 * n : n) * sizeof(T_Scalar);
            cuSolverInt sizeW = n * sizeof(T_CScalar);
            cuSolverInt sizeVL = calcLeft ? sizeA : 0;
            cuSolverInt sizeVR = calcRight ? sizeA : 0;

            void *W  = customWork().data();
            void *VL = calcLeft ? W + sizeW : nullptr;
            void *VR = calcRight ? W + sizeW + sizeVL : nullptr;
            void *vA = W + sizeW + sizeVL + sizeVR; // reuse space after W, VL, VR

            memCopyD2D<T_Scalar>(n, n, 
                                 A.values(), A.ld(), 
                                 static_cast<T_Scalar*>(vA), n);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgeev(handle(),
                            params,
                            calcLeft ? cusolverEigMode_t::CUSOLVER_EIG_MODE_VECTOR : cusolverEigMode_t::CUSOLVER_EIG_MODE_NOVECTOR,
                            calcRight ? cusolverEigMode_t::CUSOLVER_EIG_MODE_VECTOR : cusolverEigMode_t::CUSOLVER_EIG_MODE_NOVECTOR,
                            n,
                            TypeTraits<T_Scalar>::cuda_type(), vA, n, // vA is overwritten 
                            TypeTraits<T_CScalar>::cuda_type(), W,
                            TypeTraits<T_Scalar>::cuda_type(), VL, n,
                            TypeTraits<T_Scalar>::cuda_type(), VR, n,
                            TypeTraits<T_Scalar>::cuda_type(),
                            deviceWork().data(), m_workspaceInBytesOnDevice,
                            hostWork().data(), m_workspaceInBytesOnHost,
                            info().data());

            err::check_cusolver(cusolverStatus);

            // This is not needed if conjugate eigenvalues has positive imaginary part first
            // blk::dns::geev_order_eigs(n, static_cast<T_CScalar*>(W));

            m_problemCudaType = TypeTraits<T_Scalar>::cuda_type();
            m_problemDim = n;
        }

        /**
         * @brief Retrieves computed eigenvalues from the eigenvalue decomposition.
         * @details Extracts the eigenvalues computed by a previous @ref executeGeev call
         *          and stores them in the provided vector. The eigenvalues are returned
         *          as complex values.
         * @tparam T_Scalar The scalar type (real or complex).
         * @param[out] eigs Vector to store the computed eigenvalues.
         */
        template <typename T_Scalar>
        void geevGetEigenvalues(dns::CxVector<T_Scalar>& eigs) const
        {
            int_t n = m_problemDim;
            const void *W = customWork().data();
            ::cla3p::Guard<dns::CxVector<T_Scalar>> internalEigs = dns::CxVector<T_Scalar>::view(n, static_cast<const T_Scalar*>(W));
            eigs = internalEigs.get();
        }

        /**
         * @brief Retrieves computed eigenvectors from the eigenvalue decomposition.
         * @details Extracts the left and/or right eigenvectors computed by a previous
         *          @ref executeGeev call. For real matrices with complex eigenvalues,
         *          this method converts the compact real representation to full complex eigenvectors.
         * @tparam T_Scalar The scalar type (real or complex).
         * @param[in] calcLeft If true, retrieves left eigenvectors.
         * @param[in] calcRight If true, retrieves right eigenvectors.
         * @param[out] leftEigenvectors Matrix to store the left eigenvectors (if calcLeft is true).
         * @param[out] rightEigenvectors Matrix to store the right eigenvectors (if calcRight is true).
         */
        template <typename T_Scalar>
        void geevGetEigenvectors(bool calcLeft, 
                                 bool calcRight, 
                                 dns::CxMatrix<T_Scalar>& leftEigenvectors, 
                                 dns::CxMatrix<T_Scalar>& rightEigenvectors) const
        {
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            int_t n = m_problemDim;

            if(calcLeft) {
                if(!leftEigenvectors) leftEigenvectors = dns::CxMatrix<T_Scalar>(n, n);
                ::cla3p::similarity_dim_check(leftEigenvectors.nrows(), n);
                ::cla3p::similarity_dim_check(leftEigenvectors.ncols(), n);
            }

            if(calcRight) {
                if(!rightEigenvectors) rightEigenvectors = dns::CxMatrix<T_Scalar>(n, n);
                ::cla3p::similarity_dim_check(rightEigenvectors.nrows(), n);
                ::cla3p::similarity_dim_check(rightEigenvectors.ncols(), n);
            }

            cuSolverInt sizeW = n * sizeof(T_Scalar);
            cuSolverInt sizeVL = n * n * sizeof(T_RScalar);

            const void *W  = customWork().data();
            const void *VL = calcLeft ? W + sizeW : nullptr;
            const void *VR = calcRight ? W + sizeW + sizeVL : nullptr;

            if(calcLeft) {
                blk::dns::geevCalculateComplexEigenvectors(n, 
                                                           static_cast<const T_Scalar*>(W), 
                                                           static_cast<const T_RScalar*>(VL), n,
                                                           leftEigenvectors.values(), leftEigenvectors.ld());
            } // calcLeft

            if(calcRight) {
                blk::dns::geevCalculateComplexEigenvectors(n, 
                                                           static_cast<const T_Scalar*>(W), 
                                                           static_cast<const T_RScalar*>(VR), n,
                                                           rightEigenvectors.values(), rightEigenvectors.ld());
            } // calcRight
        }

    private:
        DeviceBuffer<cuSolverInt>& ipiv() { return m_ipiv; }
        DeviceBuffer<infoInt>& info() { return m_info; }
        DeviceBufferVoid& customWork() { return m_customBuffer; }
        DeviceBufferVoid& deviceWork() { return m_deviceBuffer; }
        DeviceBufferVoid& hostWork() { return m_hostBuffer; }

        const DeviceBuffer<cuSolverInt>& ipiv() const { return m_ipiv; }
        const DeviceBuffer<infoInt>& info() const { return m_info; }
        const DeviceBufferVoid& customWork() const { return m_customBuffer; }
        const DeviceBufferVoid& deviceWork() const { return m_deviceBuffer; }
        const DeviceBufferVoid& hostWork() const { return m_hostBuffer; }

    private:
        cusolverDnHandle_t m_handle{nullptr};

        size_t m_workspaceInBytesOnDevice;
        size_t m_workspaceInBytesOnHost;
        cudaDataType m_problemCudaType;
        cuSolverInt m_problemDim;

        DeviceBuffer<cuSolverInt> m_ipiv;
        DeviceBuffer<infoInt> m_info;
        DeviceBufferVoid m_customBuffer;
        DeviceBufferVoid m_deviceBuffer;
        DeviceBufferVoid m_hostBuffer;

        void defaults();
};

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_cusolver_utils
 * @brief Returns the global cuSOLVER handler instance.
 * @details Provides access to a singleton cuSOLVER handler that can be used throughout
 *          the application for cuSOLVER operations.
 * @return Reference to the global @ref CuSolverHandler instance.
 */
inline CuSolverHandler& globalCuSolverHandler()
{
    static CuSolverHandler handler;
    return handler;
}

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUSOLVER_HANDLER_HPP_
