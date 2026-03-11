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

#include <culite/types/scalar.hpp>
#include <culite/error/cuda.hpp>
#include <culite/error/exceptions.hpp>
#include <culite/support/imalloc.hpp>
#include <culite/support/utils.hpp>
#include <culite/bulk/dns1d.hpp>
#include <culite/bulk/dns2d.hpp>
#include <culite/proxies/cublas_proxy.hpp>
#include <culite/proxies/cusolver_proxy.hpp>
#include <culite/dense/dns_cxvector.hpp>
#include <culite/dense/dns_cxmatrix.hpp>

// forwards
#include <culite/checks/cla3p_forwards.hpp>

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

        class InputMeta {

            public:
                InputMeta() { defaults(); }
                InputMeta(cuSolverInt nrows, cuSolverInt ncols, cudaDataType_t cuda_type)
                    : m_nrows(nrows), m_ncols(ncols), m_cuda_type(cuda_type) {}
                
                cuSolverInt nrows() const { return m_nrows; }
                cuSolverInt ncols() const { return m_ncols; }
                cudaDataType_t cuda_type() const { return m_cuda_type; }
                
                void clear() { defaults(); }

            private:
                    cuSolverInt m_nrows;
                    cuSolverInt m_ncols;
                    cudaDataType_t m_cuda_type;

                void defaults() 
                {
                    m_nrows = 0;
                    m_ncols = 0;
                    m_cuda_type = CUDA_R_32F; // default to single precision real
                }
        };

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

            inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
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
            ::cla3p::similarity_dim_check(inputMeta().nrows(), n);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgetrs(handle(),
                             params(),
                             CUBLAS_OP_N,
                             n, B.ncols(),
                             inputMeta().cuda_type(), customWork().data(), n,
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

            inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
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

            int_t n = inputMeta().nrows();

            if(!eigs) eigs = T_Vector(n);
            ::cla3p::similarity_dim_check(eigs.size(), n);

            if(inputMeta().cuda_type() == cudaDataType::CUDA_R_32F || inputMeta().cuda_type() == cudaDataType::CUDA_R_64F) {

                using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

                T_RScalar *W = nullptr;
                geevAssignInternalPointers<T_RScalar>(false, false, n, &W, nullptr, nullptr, nullptr);

                blk::dns::geevCalculateComplexEigenvalues(n, W, eigs.values());

            } else if(inputMeta().cuda_type() == cudaDataType::CUDA_C_32F || inputMeta().cuda_type() == cudaDataType::CUDA_C_64F) {

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

            int_t n = inputMeta().nrows();

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

            if(inputMeta().cuda_type() == cudaDataType::CUDA_R_32F || inputMeta().cuda_type() == cudaDataType::CUDA_R_64F) {

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

            } else if(inputMeta().cuda_type() == cudaDataType::CUDA_C_32F || inputMeta().cuda_type() == cudaDataType::CUDA_C_64F) {

                T_Scalar *VL = nullptr;
                T_Scalar *VR = nullptr;
                geevAssignInternalPointers<T_Scalar>(calcLeft, calcRight, n, nullptr, nullptr, &VL, &VR);

                if(calcLeft)  memCopyD2D<T_Scalar>(n, n, VL, n, leftEigenvectors.values(), leftEigenvectors.ld());
                if(calcRight) memCopyD2D<T_Scalar>(n, n, VR, n, rightEigenvectors.values(), rightEigenvectors.ld());

            } else {

                throw err::CudaException("Unsupported data type for eigenvector retrieval.");

            } // problem cuda type
        }

        /**
         * @brief Reserves workspace memory for symmetric/Hermitian eigenvalue decomposition.
         * @details Computes the required workspace size and allocates buffers for performing
         *          eigenvalue decomposition on a symmetric or Hermitian matrix @p A. This includes
         *          memory for eigenvalues (real), eigenvectors (if requested), and device/host workspaces.
         * 
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property).
         * @param[in] A The symmetric/Hermitian matrix for which to reserve workspace.
         * @param[in] calcVectors If true, reserves memory for eigenvectors; if false, only eigenvalues.
         * 
         * @note This must be called before @ref executeSyevd with the same calcVectors parameter.
         * @note The matrix uplo property (upper/lower triangular storage) is preserved in workspace allocation.
         */
        template <typename T_Matrix>
        void reserveSyevd(const T_Matrix& A, bool calcVectors)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            cuSolverInt n = A.nrows();
            cuSolverInt sizeW = n * sizeof(T_RScalar);
            cuSolverInt sizeV = n * n * sizeof(T_Scalar);

            customWork().reserve(sizeW + sizeV);

            T_RScalar *W = nullptr;
            T_Scalar *V = nullptr;
            syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXsyevd_bufferSize(handle(),
                                        params(),
                                        cusolver::bool2cusolverEigMode(calcVectors),
                                        cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                                        n,
                                        TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                        TypeTraits<T_RScalar>::cuda_type(), W,
                                        TypeTraits<T_Scalar>::cuda_type(),
                                        &m_workspaceInBytesOnDevice,
                                        &m_workspaceInBytesOnHost);

            err::check_cusolver(cusolverStatus);

            info().reserve(1);
            deviceWork().reserve(m_workspaceInBytesOnDevice);
            hostWork().reserve(m_workspaceInBytesOnHost);
        }

        /**
         * @brief Computes eigenvalues and optionally eigenvectors of a symmetric/Hermitian matrix.
         * @details Performs eigenvalue decomposition on a symmetric or Hermitian matrix @p A,
         *          computing real eigenvalues @f$ \lambda @f$ and optionally eigenvectors @f$ v @f$
         *          such that @f$ A v = \lambda v @f$.
         * 
         *          This method uses cusolverDnXsyevd() internally, which employs a divide-and-conquer
         *          algorithm for efficient computation. The input matrix is copied to internal workspace
         *          before decomposition (original matrix is not modified). Results are stored internally
         *          and can be retrieved using appropriate getter methods.
         * 
         *          For symmetric/Hermitian matrices, all eigenvalues are guaranteed to be real,
         *          and eigenvectors form an orthonormal basis.
         * 
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property).
         * @param[in] A The symmetric/Hermitian matrix for which to compute eigenvalues/eigenvectors.
         * @param[in] calcVectors If true, computes eigenvectors; if false, only eigenvalues.
         * 
         * @note @ref reserveSyevd must be called first with matching calcVectors parameter.
         * @note The input matrix A is not modified; an internal copy is made.
         * @note Only the upper or lower triangular part of A is accessed, as specified by A.prop().uplo().
         */
        template <typename T_Matrix>
        void executeSyevd(const T_Matrix& A, bool calcVectors)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            cuSolverInt n = A.nrows();

            T_RScalar *W = nullptr;
            T_Scalar *V = nullptr;
            syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

            memCopyD2D<T_Scalar>(n, n, A.values(), A.ld(), V, n);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXsyevd(handle(),
                             params(),
                             cusolver::bool2cusolverEigMode(calcVectors),
                             cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                             n,
                             TypeTraits<T_Scalar>::cuda_type(), V, n,
                             TypeTraits<T_RScalar>::cuda_type(), W,
                             TypeTraits<T_Scalar>::cuda_type(),
                             deviceWork().data(), m_workspaceInBytesOnDevice,
                             hostWork().data(), m_workspaceInBytesOnHost,
                             info().data());

            err::check_cusolver(cusolverStatus);

            inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
        }

        /**
         * @brief Reserves workspace memory for selective symmetric/Hermitian eigenvalue decomposition.
         * @details Computes the required workspace size and allocates buffers for performing
         *          selective eigenvalue decomposition on a symmetric or Hermitian matrix @p A.
         *          This is an expert routine that allows computing only a subset of eigenvalues
         *          and eigenvectors based on value range or index range criteria.
         * 
         *          **Eigenvalue Range Options:**
         * 
         *          | Value | Description | Parameters Used |
         *          |-------|-------------|------------------|
         *          | All   | Computes all eigenvalues (equivalent to syevd) | - |
         *          | Value | Computes eigenvalues in the half-open interval (vl, vu] | vl, vu |
         *          | Index | Computes eigenvalues with indices il through iu (1-based indexing) | il, iu |
         * 
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property).
         * @param[in] A The symmetric/Hermitian matrix for which to reserve workspace.
         * @param[in] calcVectors If true, reserves memory for eigenvectors; if false, only eigenvalues.
         * @param[in] range Specifies the range of eigenvalues to compute (see table above).
         * @param[in] il Lower index of the eigenvalue range (1-based, used when range = Index).
         * @param[in] iu Upper index of the eigenvalue range (1-based, used when range = Index).
         * @param[in] vl Lower bound of the eigenvalue interval (used when range = Value).
         * @param[in] vu Upper bound of the eigenvalue interval (used when range = Value).
         * 
         * @note This must be called before @ref executeSyevdx with matching parameters.
         */
        template <typename T_Matrix>
        void reserveSyevdx(const T_Matrix& A, bool calcVectors, 
                           eigRange_t range, 
                           cuSolverInt il, 
                           cuSolverInt iu, 
                           typename TypeTraits<typename T_Matrix::value_type>::real_type vl, 
                           typename TypeTraits<typename T_Matrix::value_type>::real_type vu)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            cuSolverInt n = A.nrows();
            cuSolverInt sizeW = n * sizeof(T_RScalar);
            cuSolverInt sizeV = n * n * sizeof(T_Scalar);

            customWork().reserve(sizeW + sizeV);

            T_RScalar *W = nullptr;
            T_Scalar *V = nullptr;
            syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

            cuSolverInt eigsFound = 0;

            cusolverStatus_t cusolverStatus = 
            cusolverDnXsyevdx_bufferSize(handle(),
                                         params(),
                                         cusolver::bool2cusolverEigMode(calcVectors),
                                         cusolver::cla3pEigRange2cusolverEigRange(range),
                                         cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                                         n,
                                         TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                         &vl, &vu, il, iu, &eigsFound,
                                         TypeTraits<T_RScalar>::cuda_type(), W,
                                         TypeTraits<T_Scalar>::cuda_type(),
                                         &m_workspaceInBytesOnDevice,
                                         &m_workspaceInBytesOnHost);

            err::check_cusolver(cusolverStatus);

            info().reserve(1);
            deviceWork().reserve(m_workspaceInBytesOnDevice);
            hostWork().reserve(m_workspaceInBytesOnHost);
        }

        /**
         * @brief Computes selected eigenvalues and optionally eigenvectors of a symmetric/Hermitian matrix.
         * @details Performs selective eigenvalue decomposition on a symmetric or Hermitian matrix @p A,
         *          computing a subset of real eigenvalues @f$ \lambda @f$ and optionally eigenvectors @f$ v @f$
         *          such that @f$ A v = \lambda v @f$. The subset is determined by the range parameter.
         * 
         *          This method uses cusolverDnXsyevdx() internally, which employs a divide-and-conquer
         *          algorithm optimized for selective computation. The input matrix is copied to internal
         *          workspace before decomposition (original matrix is not modified). Results are stored
         *          internally and can be retrieved using appropriate getter methods.
         * 
         *          For symmetric/Hermitian matrices, all eigenvalues are guaranteed to be real,
         *          and eigenvectors form an orthonormal basis. This expert routine is more efficient
         *          than computing all eigenvalues when only a subset is needed.
         * 
         *          **Eigenvalue Range Options:**
         * 
         *          | Value | Description | Parameters Used |
         *          |-------|-------------|------------------|
         *          | All   | Behaves identically to executeSyevd and returns n eigenvalues | - |
         *          | Value | Computes eigenvalues in the half-open interval (vl, vu] | vl, vu |
         *          | Index | Eigenvalues are sorted in ascending order (1-based indexing) | il, iu |
         * 
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property).
         * @param[in] A The symmetric/Hermitian matrix for which to compute eigenvalues/eigenvectors.
         * @param[in] calcVectors If true, computes eigenvectors; if false, only eigenvalues.
         * @param[in] range Specifies the range of eigenvalues to compute (see table above).
         * @param[in] il Lower index of the eigenvalue range (1-based, used when range = Index).
         * @param[in] iu Upper index of the eigenvalue range (1-based, used when range = Index).
         * @param[in] vl Lower bound of the eigenvalue interval (used when range = Value).
         * @param[in] vu Upper bound of the eigenvalue interval (used when range = Value).
         * @return The number of eigenvalues found that satisfy the range criteria.
         * 
         * @note @ref reserveSyevdx must be called first with matching parameters.
         * @note The input matrix A is not modified; an internal copy is made.
         * @note Only the upper or lower triangular part of A is accessed, as specified by A.prop().uplo().
         */
        template <typename T_Matrix>
        int_t executeSyevdx(const T_Matrix& A, bool calcVectors, 
                            eigRange_t range, 
                            cuSolverInt il, 
                            cuSolverInt iu, 
                            typename TypeTraits<typename T_Matrix::value_type>::real_type vl, 
                            typename TypeTraits<typename T_Matrix::value_type>::real_type vu)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            cuSolverInt n = A.nrows();

            T_RScalar *W = nullptr;
            T_Scalar *V = nullptr;
            syevdAssignInternalPointers<T_Scalar>(n, &W, &V);

            memCopyD2D<T_Scalar>(n, n, A.values(), A.ld(), V, n);

            cuSolverInt eigsFound = 0;

            cusolverStatus_t cusolverStatus = 
            cusolverDnXsyevdx(handle(),
                              params(),
                              cusolver::bool2cusolverEigMode(calcVectors),
                              cusolver::cla3pEigRange2cusolverEigRange(range),
                              cublas::cla3pUplo2cublasUplo(A.prop().uplo()),
                              n,
                              TypeTraits<T_Scalar>::cuda_type(), V, n,
                              &vl, &vu, il, iu, &eigsFound,
                              TypeTraits<T_RScalar>::cuda_type(), W,
                              TypeTraits<T_Scalar>::cuda_type(),
                              deviceWork().data(), m_workspaceInBytesOnDevice,
                              hostWork().data(), m_workspaceInBytesOnHost,
                              info().data());

            err::check_cusolver(cusolverStatus);

            inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());

            return eigsFound;
        }

        /**
         * @brief Reserves workspace memory for Singular Value Decomposition (SVD).
         * @details Computes the required workspace size and allocates buffers for performing
         *          SVD on matrix @p A. This includes memory for singular values, singular vectors
         *          (U and V^T matrices), working copy of the input matrix, and device/host workspaces.
         * 
         *          **Left Singular Vectors (U) Policy:**
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Compute all m columns of U | m×m |
         *          | Limited | Compute first min(m,n) columns of U (economy-size) | m×min(m,n) |
         *          | NoCalculation | Do not compute U | - |
         * 
         *          **Right Singular Vectors (V^T) Policy:**
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Compute all n rows of V^T | n×n |
         *          | Limited | Compute first min(m,n) rows of V^T (economy-size) | min(m,n)×n |
         *          | NoCalculation | Do not compute V^T | - |
         * 
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix for which to reserve workspace (m×n).
         * @param[in] policyU Policy for computing left singular vectors U (see table above).
         * @param[in] policyVT Policy for computing right singular vectors V^T (see table above).
         * 
         * @note This must be called before @ref executeGesvd with matching policyU and policyVT parameters.
         * @note For economy-size SVD (Limited policy), memory usage is significantly reduced compared to Full.
         */
        template <typename T_Matrix>
        void reserveGesvd(const T_Matrix& A, svdPolicy_t policyU, svdPolicy_t policyVT)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            cuSolverInt m = A.nrows();
            cuSolverInt n = A.ncols();
            cuSolverInt k = std::min(m, n);

            // TODO: If allocate with Real size, cuSOLVER returns 
            //       CUDA_ERROR_MISALIGNED_ADDRESS error (for complex types only)
            // cuSolverInt sizeS = k * sizeof(T_RScalar);
            cuSolverInt sizeS = k * sizeof(T_Scalar);
            cuSolverInt sizeA = m * n * sizeof(T_Scalar);
            cuSolverInt sizeU = svdVectorSize(policyU, m, k) * sizeof(T_Scalar);
            cuSolverInt sizeV = svdVectorSize(policyVT, n, k) * sizeof(T_Scalar);

            customWork().reserve(sizeS + sizeA + sizeU + sizeV);

            T_RScalar *S = nullptr;
            T_Scalar *U = nullptr;
            T_Scalar *VT = nullptr;
            gesvdAssignInternalPointers<T_Scalar>(policyU, policyVT, m, n, &S, nullptr, &U, &VT);

            cuSolverInt ldu = m;
            cuSolverInt ldvt = (policyVT == svdPolicy_t::Limited ? k : n);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgesvd_bufferSize(handle(),
                                        params(),
                                        static_cast<signed char>(policyU),
                                        static_cast<signed char>(policyVT),
                                        m, n,
                                        TypeTraits<T_Scalar>::cuda_type(), A.values(), A.ld(),
                                        TypeTraits<T_RScalar>::cuda_type(), S,
                                        TypeTraits<T_Scalar>::cuda_type(), U, ldu,
                                        TypeTraits<T_Scalar>::cuda_type(), VT, ldvt,
                                        TypeTraits<T_Scalar>::cuda_type(),
                                        &m_workspaceInBytesOnDevice,
                                        &m_workspaceInBytesOnHost);

            err::check_cusolver(cusolverStatus);

            info().reserve(1);
            deviceWork().reserve(m_workspaceInBytesOnDevice);
            hostWork().reserve(m_workspaceInBytesOnHost);
        }

        /**
         * @brief Computes the Singular Value Decomposition (SVD) of a matrix.
         * @details Performs SVD on matrix @p A, computing the factorization @f$ A = U \Sigma V^T @f$
         *          where @f$ U @f$ and @f$ V @f$ are orthogonal/unitary matrices, and @f$ \Sigma @f$
         *          is a diagonal matrix containing the singular values in descending order.
         * 
         *          This method uses cusolverDnXgesvd() internally. The input matrix is copied to
         *          internal workspace before decomposition (original matrix is not modified).
         *          Results are stored internally and can be retrieved using appropriate getter methods.
         * 
         *          The singular values are always real and non-negative, even for complex input matrices.
         *          The computation policies control whether full or economy-size decompositions are performed.
         * 
         *          **Left Singular Vectors (U) Policy:**
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Compute all m columns of U | m×m |
         *          | Limited | Compute first min(m,n) columns of U (economy-size) | m×min(m,n) |
         *          | NoCalculation | Do not compute U | - |
         * 
         *          **Right Singular Vectors (V^T) Policy:**
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Compute all n rows of V^T | n×n |
         *          | Limited | Compute first min(m,n) rows of V^T (economy-size) | min(m,n)×n |
         *          | NoCalculation | Do not compute V^T | - |
         * 
         * @tparam T_Matrix The matrix type.
         * @param[in] A The matrix to decompose (m×n).
         * @param[in] policyU Policy for computing left singular vectors U (see table above).
         * @param[in] policyVT Policy for computing right singular vectors V^T (see table above).
         * 
         * @note @ref reserveGesvd must be called first with matching policyU and policyVT parameters.
         * @note The input matrix A is not modified; an internal copy is made.
         * @note Singular values are returned in descending order: σ₁ ≥ σ₂ ≥ ... ≥ σₘᵢₙ₍ₘ,ₙ₎ ≥ 0.
         * @note For economy-size decomposition (Limited policy), computational cost is reduced.
         */
        template <typename T_Matrix>
        void executeGesvd(const T_Matrix& A, svdPolicy_t policyU, svdPolicy_t policyVT)
        {
            using T_Scalar = typename T_Matrix::value_type;
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            cuSolverInt m = A.nrows();
            cuSolverInt n = A.ncols();
            cuSolverInt k = std::min(m, n);

            T_RScalar *S = nullptr;
            T_Scalar *vA = nullptr;
            T_Scalar *U = nullptr;
            T_Scalar *VT = nullptr;
            gesvdAssignInternalPointers<T_Scalar>(policyU, policyVT, m, n, &S, &vA, &U, &VT);

            cuSolverInt lda = m;
            cuSolverInt ldu = m;
            cuSolverInt ldvt = (policyVT == svdPolicy_t::Limited ? k : n);

            memCopyD2D<T_Scalar>(m, n, A.values(), A.ld(), vA, lda);

            cusolverStatus_t cusolverStatus = 
            cusolverDnXgesvd(handle(),
                             params(),
                             static_cast<signed char>(policyU),
                             static_cast<signed char>(policyVT),
                             m, n,
                             TypeTraits<T_Scalar>::cuda_type(), vA, lda,
                             TypeTraits<T_RScalar>::cuda_type(), S,
                             TypeTraits<T_Scalar>::cuda_type(), U, ldu,
                             TypeTraits<T_Scalar>::cuda_type(), VT, ldvt,
                             TypeTraits<T_Scalar>::cuda_type(),
                             deviceWork().data(), m_workspaceInBytesOnDevice,
                             hostWork().data(), m_workspaceInBytesOnHost,
                             info().data());

            err::check_cusolver(cusolverStatus);

            inputMeta() = InputMeta(A.nrows(), A.ncols(), TypeTraits<T_Scalar>::cuda_type());
        }

        /**
         * @brief Retrieves computed singular values from the Singular Value Decomposition.
         * @details Extracts the singular values computed by a previous @ref executeGesvd call
         *          and stores them in the provided vector. Singular values are always real and
         *          non-negative, returned in descending order: σ₁ ≥ σ₂ ≥ ... ≥ σₘᵢₙ₍ₘ,ₙ₎ ≥ 0.
         *          
         *          For an m×n matrix, the number of singular values is min(m,n).
         * @tparam T_Vector The vector type (must be a real vector type).
         * @param[out] sigma Real vector to store the computed singular values.
         * @throws CudaException if T_Vector contains complex scalars instead of real.
         * @note The output vector must have real value type (e.g., RdVector for double precision).
         * @note The vector will be automatically sized to min(m,n) if not already allocated.
         */
        template <typename T_Vector>
        void gesvdGetSingularValues(T_Vector& sigma)
        {
            using T_Scalar = typename T_Vector::value_type;

            if(!TypeTraits<T_Scalar>::is_real()) {
                throw err::CudaException("Gesvd only returns real singular values.");
            }

            int_t m = inputMeta().nrows();
            int_t n = inputMeta().ncols();
            int_t k = std::min(m, n);

            if(!sigma) sigma = T_Vector(k);
            ::cla3p::similarity_dim_check(sigma.size(), k);

            T_Scalar *S = nullptr;
            gesvdAssignInternalPointers<T_Scalar>(svdPolicy_t::NoCalculation, 
                                                  svdPolicy_t::NoCalculation, 
                                                  m, n, &S, nullptr, nullptr, nullptr);

            memCopyD2D<T_Scalar>(k, S, sigma.values());
        }

        /**
         * @brief Retrieves computed singular vectors from the Singular Value Decomposition.
         * @details Extracts the left singular vectors (U) and/or right singular vectors
         *          computed by a previous @ref executeGesvd call. The matrices returned depend
         *          on the policies specified, which must match those used in the executeGesvd call.
         * 
         *          For an m×n matrix decomposition @f$ A = U \Sigma V^T @f$, the singular vectors
         *          form orthogonal/unitary matrices satisfying specific size constraints based on
         *          the computation policy.
         * 
         *          The @p transposeVT parameter controls whether the right singular vectors are
         *          returned as @f$ V^T @f$ (transpose format, as computed by cuSOLVER) or as
         *          @f$ V @f$ (standard format, transposed for convenience).
         * 
         *          **Left Singular Vectors (U) Policy:**
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Returns all m columns of U | m×m |
         *          | Limited | Returns first min(m,n) columns of U (economy-size) | m×min(m,n) |
         *          | NoCalculation | Does not retrieve U (U is not modified) | - |
         * 
         *          **Right Singular Vectors (V/V^T) Policy:**
         * 
         *          When transposeVT = true (default), returns V:
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Returns all n columns of V | n×n |
         *          | Limited | Returns first min(m,n) columns of V (economy-size) | n×min(m,n) |
         *          | NoCalculation | Does not retrieve V (V is not modified) | - |
         * 
         *          When transposeVT = false, returns V^T:
         * 
         *          | Policy | Description | Matrix Size |
         *          |--------|-------------|-------------|
         *          | Full   | Returns all n rows of V^T | n×n |
         *          | Limited | Returns first min(m,n) rows of V^T (economy-size) | min(m,n)×n |
         *          | NoCalculation | Does not retrieve V^T (V is not modified) | - |
         * 
         * @tparam T_Matrix The matrix type.
         * @param[in] policyU Policy for retrieving left singular vectors (must match executeGesvd).
         * @param[in] policyVT Policy for retrieving right singular vectors (must match executeGesvd).
         * @param[out] U Matrix to store the left singular vectors (if policyU ≠ NoCalculation).
         * @param[out] V Matrix to store the right singular vectors, either as V or V^T depending on transposeVT.
         * @param[in] transposeVT If true (default), returns V by transposing V^T; if false, returns V^T directly.
         * 
         * @note The policyU and policyVT parameters must match those used in the executeGesvd call.
         * @note Output matrices will be automatically sized if not already allocated.
         * @note If a policy is NoCalculation, the corresponding matrix parameter is ignored.
         * @note Setting transposeVT = true requires additional transpose operation but provides V in standard format.
         */
        template <typename T_Matrix>
        void gesvdGetSingularVectors(svdPolicy_t policyU, 
                                     svdPolicy_t policyVT, 
                                     T_Matrix& U, 
                                     T_Matrix& V,
                                     bool transposeVT = true)
        {
            using T_Scalar = typename T_Matrix::value_type;

            int_t m = inputMeta().nrows();
            int_t n = inputMeta().ncols();
            int_t k = std::min(m, n);

            int_t nrowsU = 0;
            int_t ncolsU = 0;
            int_t nrowsVT = 0;
            int_t ncolsVT = 0;

            if(policyU != svdPolicy_t::NoCalculation) {
                nrowsU = m;
                ncolsU = (policyU == svdPolicy_t::Limited ? k : m);
                if(!U) U = T_Matrix(nrowsU, ncolsU);
                ::cla3p::similarity_dim_check(U.nrows(), nrowsU);
                ::cla3p::similarity_dim_check(U.ncols(), ncolsU);
            }

            if(policyVT != svdPolicy_t::NoCalculation) {
                nrowsVT = (policyVT == svdPolicy_t::Limited ? k : n);
                ncolsVT = n;
                int_t nrowsV = transposeVT ? ncolsVT : nrowsVT;
                int_t ncolsV = transposeVT ? nrowsVT : ncolsVT;
                if(!V) V = T_Matrix(nrowsV, ncolsV);
                ::cla3p::similarity_dim_check(V.nrows(), nrowsV);
                ::cla3p::similarity_dim_check(V.ncols(), ncolsV);
            }

            T_Scalar *storedU = nullptr;
            T_Scalar *storedVT = nullptr;
            gesvdAssignInternalPointers<T_Scalar>(policyU, policyVT, m, n, nullptr, nullptr, &storedU, &storedVT);

            if(policyU != svdPolicy_t::NoCalculation) {
                memCopyD2D<T_Scalar>(nrowsU, ncolsU, storedU, nrowsU, U.values(), U.ld());
            } // policyU
            
            if(policyVT != svdPolicy_t::NoCalculation) {
                if(transposeVT) {
                    blk::dns::ctranspose<T_Scalar>(nrowsVT, ncolsVT, storedVT, nrowsVT, V.values(), V.ld());
                } else {
                    memCopyD2D<T_Scalar>(nrowsVT, ncolsVT, storedVT, nrowsVT, V.values(), V.ld());
                }
            } // policyVT
        }




        // ------------ Helper Methods for Internal Pointer Management -----------

        cuSolverInt svdVectorSize(svdPolicy_t policy, cuSolverInt n, cuSolverInt k) const
        {
            switch(policy) {
                case svdPolicy_t::Full         : return n * n;
                case svdPolicy_t::Limited      : return n * k;
                case svdPolicy_t::NoCalculation: return 0;
                default: return 0;
            }
        }

    private:
        InputMeta& inputMeta() { return m_inputMeta; }
        DeviceBuffer<cuSolverInt>& ipiv() { return m_ipiv; }
        DeviceBuffer<infoInt>& info() { return m_info; }
        DeviceBufferVoid& customWork() { return m_customBuffer; }
        DeviceBufferVoid& deviceWork() { return m_deviceBuffer; }
        PinnedBufferVoid& hostWork() { return m_hostBuffer; }

        const InputMeta& inputMeta() const { return m_inputMeta; }
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

        template <typename T_Scalar>
        void syevdAssignInternalPointers(int_t n, typename TypeTraits<T_Scalar>::real_type** W, T_Scalar** V)
        {
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            if(W) *W = nullptr;
            if(V) *V = nullptr;

            char *charBuffer = static_cast<char*>(customWork().data());

            std::size_t sizeW = n * sizeof(T_RScalar);
            std::size_t sizeV = n * n * sizeof(T_Scalar);

            if(W) { *W = reinterpret_cast<T_RScalar*>(charBuffer); } charBuffer += sizeW;
            if(V) { *V = reinterpret_cast<T_Scalar*>(charBuffer); } charBuffer += sizeV;
        }

        template <typename T_Scalar>
        void gesvdAssignInternalPointers(svdPolicy_t policyU, 
                                         svdPolicy_t policyV,
                                         int_t m,
                                         int_t n, 
                                         typename TypeTraits<T_Scalar>::real_type** S, 
                                         T_Scalar** A,
                                         T_Scalar** U, 
                                         T_Scalar** V)
        {
            using T_RScalar = typename TypeTraits<T_Scalar>::real_type;

            int_t k = std::min(m, n);

            if(S) *S = nullptr;
            if(A) *A = nullptr;
            if(U) *U = nullptr;
            if(V) *V = nullptr;

            char *charBuffer = static_cast<char*>(customWork().data());

            // TODO: If allocate with Real size, cuSOLVER returns 
            //       CUDA_ERROR_MISALIGNED_ADDRESS error (for complex types only)
            // cuSolverInt sizeS = k * sizeof(T_RScalar);
            cuSolverInt sizeS = k * sizeof(T_Scalar);
            std::size_t sizeA = m * n * sizeof(T_Scalar);
            std::size_t sizeU = svdVectorSize(policyU, m, k) * sizeof(T_Scalar);
            std::size_t sizeV = svdVectorSize(policyV, n, k) * sizeof(T_Scalar);

            if(S) { *S = reinterpret_cast<T_RScalar*>(charBuffer); } charBuffer += sizeS;
            if(A) { *A = reinterpret_cast<T_Scalar *>(charBuffer); } charBuffer += sizeA;
            if(U) { *U = reinterpret_cast<T_Scalar *>(charBuffer); } charBuffer += sizeU;
            if(V) { *V = reinterpret_cast<T_Scalar *>(charBuffer); } charBuffer += sizeV;
        }

    private:
        cusolverDnHandle_t m_handle{nullptr};
        cusolver::DnParams m_params;

        size_t m_workspaceInBytesOnDevice;
        size_t m_workspaceInBytesOnHost;
        InputMeta m_inputMeta;
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
