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

#include <culite/types/integer.hpp>
#include <culite/types/scalar.hpp>
#include <culite/support/imalloc.hpp>
#include <culite/proxies/cusolver_proxy.hpp>

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
                InputMeta();
                InputMeta(cuSolverInt nrows, cuSolverInt ncols, cudaDataType_t cuda_type);

                cuSolverInt nrows() const;
                cuSolverInt ncols() const;
                cudaDataType_t cuda_type() const;

                void clear();

            private:
                cuSolverInt m_nrows;
                cuSolverInt m_ncols;
                cudaDataType_t m_cuda_type;

                void defaults();
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
        cusolverDnHandle_t handle();

        /**
         * @brief Get the cuSolver parameters object.
         * @return The cuSolver DN parameters handle.
         */
        cusolverDnParams_t params() const;

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
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The matrix for which to reserve workspace.
         */
        template <typename T_Matrix>
        void reserveLU(const T_Matrix& A);

        /**
         * @brief Performs LU factorization with partial pivoting.
         * @details Computes the LU decomposition @f$ P A = L U @f$ of matrix @p A, where
         *          @p P is a permutation matrix, @p L is lower triangular, and @p U is upper triangular.
         *          The factorization is stored internally for subsequent solve operations.
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The matrix to factorize.
         */
        template <typename T_Matrix>
        void decomposeLU(const T_Matrix& A);

        /**
         * @brief Solves a linear system using the computed LU factorization.
         * @details Solves the system @f$ A X = B @f$ using the previously computed LU factorization.
         *          The solution is stored in-place in matrix @p B.
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in,out] B On input, the right-hand side matrix; on output, the solution matrix.
         */
        template <typename T_Matrix>
        void solveLU(T_Matrix& B);

        /**
         * @brief Reserves workspace memory for eigenvalue decomposition.
         * @details Computes the required workspace size and allocates buffers for performing
         *          eigenvalue decomposition on matrix @p A. This includes memory for eigenvalues,
         *          eigenvectors (if requested), and device/host workspaces.
         * 
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The matrix for which to reserve workspace.
         * @param[in] calcLeft If true, reserves memory for left eigenvectors.
         * @param[in] calcRight If true, reserves memory for right eigenvectors.
         * 
         * @note This must be called before @ref executeGeev with the same calcLeft/calcRight parameters.
         */
        template <typename T_Matrix>
        void reserveGeev(const T_Matrix& A, bool calcLeft, bool calcRight);

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
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The matrix for which to compute eigenvalues/eigenvectors.
         * @param[in] calcLeft If true, computes left eigenvectors.
         * @param[in] calcRight If true, computes right eigenvectors.
         * 
         * @note @ref reserveGeev must be called first with matching parameters.
         * @note The input matrix A is not modified; an internal copy is made.
         */
        template <typename T_Matrix>
        void executeGeev(const T_Matrix& A, bool calcLeft, bool calcRight);

        /**
         * @brief Retrieves computed eigenvalues from the eigenvalue decomposition.
         * @details Extracts the eigenvalues computed by a previous @ref executeGeev call
         *          and stores them in the provided vector. Eigenvalues are always returned
         *          as complex values, even for real input matrices.
         * @tparam T_Vector The vector type. Supported: dns::CdVector, dns::CfVector.
         * @param[out] eigs Complex vector to store the computed eigenvalues.
         * @throws CudaException if T_Vector contains real scalars instead of complex.
         * @note The output vector must have complex value type (e.g., dns::CdVector for double precision).
         */
        template <typename T_Vector>
        void geevGetEigenvalues(T_Vector& eigs);

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
         * @tparam T_Matrix The matrix type. Supported: dns::CdMatrix, dns::CfMatrix.
         * @param[in] calcLeft If true, retrieves left eigenvectors (must match executeGeev parameters).
         * @param[in] calcRight If true, retrieves right eigenvectors (must match executeGeev parameters).
         * @param[out] leftEigenvectors Complex matrix to store the left eigenvectors (if calcLeft is true).
         * @param[out] rightEigenvectors Complex matrix to store the right eigenvectors (if calcRight is true).
         * 
         * @note The output matrices must have complex value type (e.g., dns::CdMatrix for double precision).
         * @note The calcLeft and calcRight parameters must match those used in the executeGeev call.
         */
        template <typename T_Matrix>
        void geevGetEigenvectors(bool calcLeft, 
                                 bool calcRight, 
                                 T_Matrix& leftEigenvectors, 
                                 T_Matrix& rightEigenvectors);

        /**
         * @brief Reserves workspace memory for symmetric/Hermitian eigenvalue decomposition.
         * @details Computes the required workspace size and allocates buffers for performing
         *          eigenvalue decomposition on a symmetric or Hermitian matrix @p A. This includes
         *          memory for eigenvalues (real), eigenvectors (if requested), and device/host workspaces.
         * 
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property). Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The symmetric/Hermitian matrix for which to reserve workspace.
         * @param[in] calcVectors If true, reserves memory for eigenvectors; if false, only eigenvalues.
         * 
         * @note This must be called before @ref executeSyevd with the same calcVectors parameter.
         * @note The matrix uplo property (upper/lower triangular storage) is preserved in workspace allocation.
         */
        template <typename T_Matrix>
        void reserveSyevd(const T_Matrix& A, bool calcVectors);

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
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property). Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The symmetric/Hermitian matrix for which to compute eigenvalues/eigenvectors.
         * @param[in] calcVectors If true, computes eigenvectors; if false, only eigenvalues.
         * 
         * @note @ref reserveSyevd must be called first with matching calcVectors parameter.
         * @note The input matrix A is not modified; an internal copy is made.
         * @note Only the upper or lower triangular part of A is accessed, as specified by A.prop().uplo().
         */
        template <typename T_Matrix>
        void executeSyevd(const T_Matrix& A, bool calcVectors);

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
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property). Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
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
                           typename TypeTraits<typename T_Matrix::value_type>::real_type vu);

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
         * @tparam T_Matrix The matrix type (must have symmetric or Hermitian property). Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
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
                            typename TypeTraits<typename T_Matrix::value_type>::real_type vu);

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
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
         * @param[in] A The matrix for which to reserve workspace (m×n).
         * @param[in] policyU Policy for computing left singular vectors U (see table above).
         * @param[in] policyVT Policy for computing right singular vectors V^T (see table above).
         * 
         * @note This must be called before @ref executeGesvd with matching policyU and policyVT parameters.
         * @note For economy-size SVD (Limited policy), memory usage is significantly reduced compared to Full.
         */
        template <typename T_Matrix>
        void reserveGesvd(const T_Matrix& A, svdPolicy_t policyU, svdPolicy_t policyVT);

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
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
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
        void executeGesvd(const T_Matrix& A, svdPolicy_t policyU, svdPolicy_t policyVT);

        /**
         * @brief Retrieves computed singular values from the Singular Value Decomposition.
         * @details Extracts the singular values computed by a previous @ref executeGesvd call
         *          and stores them in the provided vector. Singular values are always real and
         *          non-negative, returned in descending order: σ₁ ≥ σ₂ ≥ ... ≥ σₘᵢₙ₍ₘ,ₙ₎ ≥ 0.
         *          
         *          For an m×n matrix, the number of singular values is min(m,n).
         * @tparam T_Vector The vector type. Supported: dns::RdVector, dns::RfVector.
         * @param[out] sigma Real vector to store the computed singular values.
         * @throws CudaException if T_Vector contains complex scalars instead of real.
         * @note The output vector must have real value type (e.g., dns::RdVector for double precision).
         * @note The vector will be automatically sized to min(m,n) if not already allocated.
         */
        template <typename T_Vector>
        void gesvdGetSingularValues(T_Vector& sigma);

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
         * @tparam T_Matrix The matrix type. Supported: dns::RdMatrix, dns::RfMatrix, dns::CdMatrix, dns::CfMatrix.
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
                                     bool transposeVT = true);




        // ------------ Helper Methods for Internal Pointer Management -----------

        cuSolverInt svdVectorSize(svdPolicy_t policy, cuSolverInt n, cuSolverInt k) const;

    private:
        InputMeta& inputMeta();
        DeviceBuffer<cuSolverInt>& ipiv();
        DeviceBuffer<infoInt>& info();
        DeviceBufferVoid& customWork();
        DeviceBufferVoid& deviceWork();
        PinnedBufferVoid& hostWork();

        const InputMeta& inputMeta() const;
        const DeviceBuffer<cuSolverInt>& ipiv() const;
        const DeviceBuffer<infoInt>& info() const;
        const DeviceBufferVoid& customWork() const;
        const DeviceBufferVoid& deviceWork() const;
        const PinnedBufferVoid& hostWork() const;

        template <typename T_Scalar>
        void geevAssignInternalPointers(bool calcLeft, 
                                        bool calcRight, 
                                        int_t n,
                                        T_Scalar** W, 
                                        T_Scalar** A,
                                        T_Scalar** VL, 
                                        T_Scalar** VR);

        template <typename T_Scalar>
        void syevdAssignInternalPointers(int_t n, typename TypeTraits<T_Scalar>::real_type** W, T_Scalar** V);

        template <typename T_Scalar>
        void gesvdAssignInternalPointers(svdPolicy_t policyU, 
                                         svdPolicy_t policyV,
                                         int_t m,
                                         int_t n, 
                                         typename TypeTraits<T_Scalar>::real_type** S, 
                                         T_Scalar** A,
                                         T_Scalar** U, 
                                         T_Scalar** V);

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
