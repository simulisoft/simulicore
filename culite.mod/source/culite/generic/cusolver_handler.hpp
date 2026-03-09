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

#include "culite/types/scalar.hpp"
#include "culite/error/cuda.hpp"
#include "culite/error/exceptions.hpp"
#include "culite/support/imalloc.hpp"
#include "culite/support/utils.hpp"
#include "culite/bulk/dns1d.hpp"
#include "culite/bulk/dns2d.hpp"
#include "culite/proxies/cusolver_proxy.hpp"
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
 * 
 *          The handler uses both device memory buffers and pinned host memory buffers
 *          to optimize performance and meet cuSOLVER requirements.
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
         * @brief Get the cuSolver parameters object.
         * @return The cuSolver DN parameters handle.
         */
        cusolverDnParams_t params() const { return m_params.get(); }

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

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrf_bufferSize(handle(),
                                        params(),
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

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrf(handle(),
                             params(),
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
            using T_Scalar = typename T_Matrix::value_type;

            cuSolverInt n = B.nrows();
            ::cla3p::similarity_dim_check(m_problemDim, n);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrs(handle(),
                             params(),
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
         * 
         *          Memory layout in customWork buffer: W | A | VL | VR
         *          - W: eigenvalues (complex)
         *          - A: working copy of input matrix
         *          - VL: left eigenvectors (if calcLeft = true)
         *          - VR: right eigenvectors (if calcRight = true)
         * 
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix for which to reserve workspace.
         * @param[in] calcLeft If true, reserves memory for left eigenvectors.
         * @param[in] calcRight If true, reserves memory for right eigenvectors.
         * 
         * @note This must be called before @ref executeGeev with the same calcLeft/calcRight parameters.
         */
        template <typename T_Matrix>
        void reserveGeev(const T_Matrix& A, bool calcLeft, bool calcRight)
        {
            using T_Scalar = typename T_Matrix::value_type;
            bool realCase = TypeTraits<T_Scalar>::is_real();

            cuSolverInt n = A.nrows();
            cuSolverInt sizeW = (realCase ? 2 * n : n) * sizeof(T_Scalar);
            cuSolverInt sizeN2 = n * n * sizeof(T_Scalar);
            cuSolverInt sizevA = sizeN2;
            cuSolverInt sizeVL = calcLeft ? sizeN2 : 0;
            cuSolverInt sizeVR = calcRight ? sizeN2 : 0;

            customWork().reserve(sizeW + sizeVL + sizeVR + sizevA);

            T_Scalar *W  = nullptr;
            T_Scalar *VL = nullptr;
            T_Scalar *VR = nullptr;
            geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, &W, nullptr, &VL, &VR);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgeev_bufferSize(handle(),
                                       params(),
                                       cusolver::bool2cusolverEigMode(calcLeft),
                                       cusolver::bool2cusolverEigMode(calcRight),
                                       n,
                                       TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                       TypeTraits<T_Scalar>::cuda_type(), W,
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
         * 
         *          This method uses cusolverDnXgeev() internally. The input matrix is copied
         *          to internal workspace before decomposition (original matrix is not modified).
         *          Results are stored internally and can be retrieved using @ref geevGetEigenvalues
         *          and @ref geevGetEigenvectors.
         * 
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix for which to compute eigenvalues/eigenvectors.
         * @param[in] calcLeft If true, computes left eigenvectors.
         * @param[in] calcRight If true, computes right eigenvectors.
         * 
         * @note @ref reserveGeev must be called first with matching parameters.
         * @note The input matrix A is not modified; an internal copy is made.
         */
        template <typename T_Matrix>
        void executeGeev(const T_Matrix& A, bool calcLeft, bool calcRight)
        {
            using T_Scalar = typename T_Matrix::value_type;

            cuSolverInt n = A.nrows();

            T_Scalar *W  = nullptr;
            T_Scalar *vA = nullptr;
            T_Scalar *VL = nullptr;
            T_Scalar *VR = nullptr;
            geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, &W, &vA, &VL, &VR);

            memCopyD2D<T_Scalar>(n, n, A.values(), A.ld(), vA, n);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgeev(handle(),
                            params(),
                            cusolver::bool2cusolverEigMode(calcLeft),
                            cusolver::bool2cusolverEigMode(calcRight),
                            n,
                            TypeTraits<T_Scalar>::cuda_type(), vA, n, // vA is overwritten 
                            TypeTraits<T_Scalar>::cuda_type(), W,
                            TypeTraits<T_Scalar>::cuda_type(), VL, n,
                            TypeTraits<T_Scalar>::cuda_type(), VR, n,
                            TypeTraits<T_Scalar>::cuda_type(),
                            deviceWork().data(), m_workspaceInBytesOnDevice,
                            hostWork().data(), m_workspaceInBytesOnHost,
                            info().data());

            err::check_cusolver(cusolverStatus);

            m_problemCudaType = TypeTraits<T_Scalar>::cuda_type();
            m_problemDim = n;
        }

        /**
         * @brief Retrieves computed eigenvalues from the eigenvalue decomposition.
         * @details Extracts the eigenvalues computed by a previous @ref executeGeev call
         *          and stores them in the provided vector. Eigenvalues are always returned
         *          as complex values, even for real input matrices.
         * @tparam T_Vector The vector type (must be a complex vector type).
         * @param[out] eigs Complex vector to store the computed eigenvalues.
         * @throws CudaException if T_Vector contains real scalars instead of complex.
         * @note The output vector must have complex value type (e.g., CdVector for double precision).
         */
        template <typename T_Vector>
        void geevGetEigenvalues(T_Vector& eigs)
        {
            using T_Scalar = typename T_Vector::value_type;

            if(TypeTraits<T_Scalar>::is_real()) {
                throw err::CudaException("Geev only returns complex eigenpairs.");
            }

            int_t n = m_problemDim;

            if(!eigs) eigs = T_Vector(n);
            ::cla3p::similarity_dim_check(eigs.size(), n);

            if(m_problemCudaType == cudaDataType::CUDA_R_32F || m_problemCudaType == cudaDataType::CUDA_R_64F) {

                using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

                T_RScalar *W = nullptr;
                geevAssignInternalPointers<T_RScalar>(false, false, n, &W, nullptr, nullptr, nullptr);

                blk::dns::geevCalculateComplexEigenvalues(n, W, eigs.values());

            } else if(m_problemCudaType == cudaDataType::CUDA_C_32F || m_problemCudaType == cudaDataType::CUDA_C_64F) {

                T_Scalar *W = nullptr;
                geevAssignInternalPointers<T_Scalar>(false, false, n, &W, nullptr, nullptr, nullptr);

                memCopyD2D<T_Scalar>(n, W, eigs.values());

            } else {

                throw err::CudaException("Unsupported data type for eigenvalue retrieval.");

            } // problem cuda type
        }

        /**
         * @brief Retrieves computed eigenvectors from the eigenvalue decomposition.
         * @details Extracts the left and/or right eigenvectors computed by a previous
         *          @ref executeGeev call. Eigenvectors are always returned as complex matrices.
         *          
         *          For real input matrices with complex eigenvalues, this method automatically
         *          converts the compact real representation (where conjugate pairs share storage)
         *          to full complex eigenvectors. For complex input matrices, eigenvectors are
         *          returned directly without conversion.
         * 
         * @tparam T_Matrix The matrix type (must be a complex matrix type).
         * @param[in] calcLeft If true, retrieves left eigenvectors (must match executeGeev parameters).
         * @param[in] calcRight If true, retrieves right eigenvectors (must match executeGeev parameters).
         * @param[out] leftEigenvectors Complex matrix to store the left eigenvectors (if calcLeft is true).
         * @param[out] rightEigenvectors Complex matrix to store the right eigenvectors (if calcRight is true).
         * 
         * @note The output matrices must have complex value type (e.g., CdMatrix for double precision).
         * @note The calcLeft and calcRight parameters must match those used in the executeGeev call.
         */
        template <typename T_Matrix>
        void geevGetEigenvectors(bool calcLeft, 
                                 bool calcRight, 
                                 T_Matrix& leftEigenvectors, 
                                 T_Matrix& rightEigenvectors)
        {
            using T_Scalar = typename T_Matrix::value_type;

            if(TypeTraits<T_Scalar>::is_real()) {
                throw err::CudaException("Geev only returns complex eigenpairs.");
            }

            int_t n = m_problemDim;

            if(calcLeft) {
                if(!leftEigenvectors) leftEigenvectors = T_Matrix(n, n);
                ::cla3p::similarity_dim_check(leftEigenvectors.nrows(), n);
                ::cla3p::similarity_dim_check(leftEigenvectors.ncols(), n);
            }

            if(calcRight) {
                if(!rightEigenvectors) rightEigenvectors = T_Matrix(n, n);
                ::cla3p::similarity_dim_check(rightEigenvectors.nrows(), n);
                ::cla3p::similarity_dim_check(rightEigenvectors.ncols(), n);
            }

            if(m_problemCudaType == cudaDataType::CUDA_R_32F || m_problemCudaType == cudaDataType::CUDA_R_64F) {

                using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

                T_RScalar *W  = nullptr;
                T_RScalar *VL = nullptr;
                T_RScalar *VR = nullptr;
                geevAssignInternalPointers<T_RScalar>(calcLeft, calcRight, n, &W, nullptr, &VL, &VR);

                if(calcLeft) {
                    blk::dns::geevCalculateComplexEigenvectors(n, W, VL, n,
                                                               leftEigenvectors.values(), 
                                                               leftEigenvectors.ld());
                } // calcLeft

                if(calcRight) {
                    blk::dns::geevCalculateComplexEigenvectors(n, W, VR, n,
                                                               rightEigenvectors.values(), 
                                                               rightEigenvectors.ld());
                } // calcRight

            } else if(m_problemCudaType == cudaDataType::CUDA_C_32F || m_problemCudaType == cudaDataType::CUDA_C_64F) {

                T_Scalar *VL = nullptr;
                T_Scalar *VR = nullptr;
                geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, nullptr, nullptr, &VL, &VR);

                if(calcLeft)  memCopyD2D<T_Scalar>(n, n, VL, n, leftEigenvectors.values(), leftEigenvectors.ld());
                if(calcRight) memCopyD2D<T_Scalar>(n, n, VR, n, rightEigenvectors.values(), rightEigenvectors.ld());

            } else {

                throw err::CudaException("Unsupported data type for eigenvector retrieval.");

            } // problem cuda type
        }

    private:
        DeviceBuffer<cuSolverInt>& ipiv() { return m_ipiv; }
        DeviceBuffer<infoInt>& info() { return m_info; }
        DeviceBufferVoid& customWork() { return m_customBuffer; }
        DeviceBufferVoid& deviceWork() { return m_deviceBuffer; }
        PinnedBufferVoid& hostWork() { return m_hostBuffer; }

        const DeviceBuffer<cuSolverInt>& ipiv() const { return m_ipiv; }
        const DeviceBuffer<infoInt>& info() const { return m_info; }
        const DeviceBufferVoid& customWork() const { return m_customBuffer; }
        const DeviceBufferVoid& deviceWork() const { return m_deviceBuffer; }
        const PinnedBufferVoid& hostWork() const { return m_hostBuffer; }

        template <typename T_Scalar>
        void geevAssignInternalPointers(bool calcLeft, 
                                        bool calcRight, 
                                        int_t n,
                                        T_Scalar** W, 
                                        T_Scalar** A,
                                        T_Scalar** VL, 
                                        T_Scalar** VR)
        {
            bool realCase = TypeTraits<T_Scalar>::is_real();

            if(W)  *W  = nullptr;
            if(A)  *A  = nullptr;
            if(VL) *VL = nullptr;
            if(VR) *VR = nullptr;
            
            char *charBuffer = static_cast<char*>(customWork().data());

            std::size_t sizeW  = (realCase ? 2 * n : n) * sizeof(T_Scalar);
            std::size_t sizeN2 = n * n * sizeof(T_Scalar);
            std::size_t sizeA  = sizeN2;
            std::size_t sizeVL = calcLeft ? sizeN2 : 0;
            std::size_t sizeVR = calcRight ? sizeN2 : 0;

            if(W)               { *W  = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeW ;
            if(A)               { *A  = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeA ;
            if(calcLeft  && VL) { *VL = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeVL;
            if(calcRight && VR) { *VR = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeVR;
        }

    private:
        cusolverDnHandle_t m_handle{nullptr};
        cusolver::DnParams m_params;

        size_t m_workspaceInBytesOnDevice;
        size_t m_workspaceInBytesOnHost;
        cudaDataType m_problemCudaType;
        cuSolverInt m_problemDim;

        DeviceBuffer<cuSolverInt> m_ipiv;
        DeviceBuffer<infoInt> m_info;
        DeviceBufferVoid m_customBuffer;
        DeviceBufferVoid m_deviceBuffer;
        PinnedBufferVoid m_hostBuffer;

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
