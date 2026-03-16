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

#ifndef CCULITE_CUSPARSE_HANDLER_HPP_
#define CCULITE_CUSPARSE_HANDLER_HPP_

/**
 * @file
 */

#include <culite/support/imalloc.hpp>
#include <culite/proxies/cusparse_proxy.hpp>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

class CuSparseHandler;

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief RAII helper for managing cuSPARSE pointer mode.
 * @details This class provides a scoped guard that temporarily sets a cuSPARSE pointer mode
 *          and automatically restores the previous mode when the object goes out of scope.
 *          This ensures proper cleanup even in the presence of exceptions.
 */
class CuSparseSpm {
    public:

        // disable copying
        CuSparseSpm(const CuSparseSpm&) = delete;
        CuSparseSpm& operator=(const CuSparseSpm&) = delete;

        /**
         * @brief Constructor.
         * @details Sets the pointer mode for the given CuSparseHandler and stores the previous mode.
         * @param[in] cuSparseHandler The CuSparseHandler whose pointer mode will be modified.
         * @param[in] mode The new pointer mode to set.
         */
        CuSparseSpm(CuSparseHandler& cuSparseHandler, cusparsePointerMode_t mode);
        
        /**
         * @brief Destructor.
         * @details Restores the pointer mode that was active before the constructor was called.
         */
        ~CuSparseSpm();

    private:
        CuSparseHandler& m_handler;
        cusparsePointerMode_t m_oldMode;
};  

/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The cuSparse handler class.
 * @details This class provides a wrapper around the cuSPARSE library for performing
 *          sparse linear algebra operations on GPU devices. It manages the cuSPARSE handle
 *          and internal workspace buffers required for factorization and solve operations.
 */
class CuSparseHandler {

    public:

        // disable copying
        CuSparseHandler(const CuSparseHandler&) = delete;
        CuSparseHandler& operator=(const CuSparseHandler&) = delete;

        /**
         * @brief Constructor.
         * @details Initializes the cuSPARSE handle and internal state.
         */
        CuSparseHandler();

        /**
         * @brief Destructor.
         * @details Destroys the cuSPARSE handle and releases all allocated resources.
         */
        ~CuSparseHandler();
        
        /**
         * @brief Get the cuSPARSE handle.
         * @return The cuSPARSE handle.
         */
        cusparseHandle_t handle();

        /**
         * @brief Clear the internal workspace buffers.
         * @details Releases all internal workspace memory allocated for cuSPARSE operations.
         */
        void clear();

        /**
         * @brief Set the pointer mode for cuSPARSE operations.
         * @details Changes the pointer mode and returns the previous mode. The pointer mode
         *          determines whether scalar values (alpha, beta) are passed by reference
         *          on the host or device.
         * @param[in] mode The new pointer mode to set.
         * @return The pointer mode that was active before this call.
         */
        cusparsePointerMode_t setPointerMode(cusparsePointerMode_t mode);
        
        /**
         * @brief Get the current pointer mode.
         * @details Returns the current pointer mode setting for cuSPARSE operations.
         * @return The current pointer mode.
         */
        cusparsePointerMode_t pointerMode();

        /**
         * @brief Reserve workspace memory for sparse matrix-vector multiplication (SpMV).
         * @details Queries the required buffer size for the SpMV operation and allocates
         *          the necessary workspace memory on the device.
         * @tparam T_Scalar The scalar type (e.g., float, double, complex).
         * @param[in] opA Operation to apply to matrix A (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] alpha Pointer to scalar alpha.
         * @param[in] matA Sparse matrix.
         * @param[in] vecX Input dense vector X.
         * @param[in] beta Pointer to scalar beta.
         * @param[in,out] vecY Dense vector Y for buffer size computation.
         * @param[in] alg Algorithm to use for SpMV (default: CUSPARSE_SPMV_CSR_ALG1).
         */
        template <typename T_Scalar>
        void reserveSpmv(op_t                             opA,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnVec<T_Scalar>& vecX,
                         const T_Scalar*                  beta,
                         cusparse::DnVec<T_Scalar>&       vecY,
                         cusparseSpMVAlg_t                alg = cusparseSpMVAlg_t::CUSPARSE_SPMV_CSR_ALG1);

        /**
         * @brief Preprocess the sparse matrix-vector multiplication operation.
         * @details Performs any necessary preprocessing steps for the SpMV operation,
         *          which may optimize subsequent repeated operations with the same sparsity pattern.
         * @tparam T_Scalar The scalar type (e.g., float, double, complex).
         * @param[in] opA Operation to apply to matrix A (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] alpha Pointer to scalar alpha.
         * @param[in] matA Sparse matrix.
         * @param[in] vecX Input dense vector X.
         * @param[in] beta Pointer to scalar beta.
         * @param[in,out] vecY Dense vector Y for preprocessing.
         * @param[in] alg Algorithm to use for SpMV (default: CUSPARSE_SPMV_CSR_ALG1).
         */
        template <typename T_Scalar>
        void preprocessSpmv(op_t                             opA,
                            const T_Scalar*                  alpha,
                            const cusparse::SpMatBase&       matA,
                            const cusparse::DnVec<T_Scalar>& vecX,
                            const T_Scalar*                  beta,
                            cusparse::DnVec<T_Scalar>&       vecY,
                            cusparseSpMVAlg_t                alg = cusparseSpMVAlg_t::CUSPARSE_SPMV_CSR_ALG1);

        /**
         * @brief Perform sparse matrix-vector multiplication.
         * @details Computes Y = alpha * op(A) * X + beta * Y, where A is a sparse matrix,
         *          X and Y are dense vectors, and op(A) is either A, A^T, or A^H.
         * @tparam T_Scalar The scalar type (e.g., float, double, complex).
         * @param[in] opA Operation to apply to matrix A (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] alpha Pointer to scalar alpha.
         * @param[in] matA Sparse matrix.
         * @param[in] vecX Input dense vector X.
         * @param[in] beta Pointer to scalar beta.
         * @param[in,out] vecY Output dense vector Y.
         * @param[in] alg Algorithm to use for SpMV (default: CUSPARSE_SPMV_CSR_ALG1).
         */
        template <typename T_Scalar>
        void performSpmv(op_t                             opA,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnVec<T_Scalar>& vecX,
                         const T_Scalar*                  beta,
                         cusparse::DnVec<T_Scalar>&       vecY,
                         cusparseSpMVAlg_t                alg = cusparseSpMVAlg_t::CUSPARSE_SPMV_CSR_ALG1);

        /**
         * @brief Reserve workspace memory for sparse matrix-matrix multiplication (SpMM).
         * @details Queries the required buffer size for the SpMM operation and allocates
         *          the necessary workspace memory on the device.
         * @tparam T_Scalar The scalar type (e.g., float, double, complex).
         * @param[in] opA Operation to apply to matrix A (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] opB Operation to apply to matrix B (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] alpha Pointer to scalar alpha.
         * @param[in] matA Sparse matrix A.
         * @param[in] matB Dense matrix B.
         * @param[in] beta Pointer to scalar beta.
         * @param[in,out] matC Dense matrix C for buffer size computation.
         * @param[in] alg Algorithm to use for SpMM (default: CUSPARSE_SPMM_CSR_ALG1).
         */
        template <typename T_Scalar>
        void reserveSpmm(op_t                             opA,
                         op_t                             opB,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnMat<T_Scalar>& matB,
                         const T_Scalar*                  beta,
                         cusparse::DnMat<T_Scalar>&       matC,
                         cusparseSpMMAlg_t                alg = cusparseSpMMAlg_t::CUSPARSE_SPMM_CSR_ALG1);

        /**
         * @brief Preprocess the sparse matrix-matrix multiplication operation.
         * @details Performs any necessary preprocessing steps for the SpMM operation,
         *          which may optimize subsequent repeated operations with the same sparsity pattern.
         * @tparam T_Scalar The scalar type (e.g., float, double, complex).
         * @param[in] opA Operation to apply to matrix A (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] opB Operation to apply to matrix B (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] alpha Pointer to scalar alpha.
         * @param[in] matA Sparse matrix A.
         * @param[in] matB Dense matrix B.
         * @param[in] beta Pointer to scalar beta.
         * @param[in,out] matC Dense matrix C for preprocessing.
         * @param[in] alg Algorithm to use for SpMM (default: CUSPARSE_SPMM_CSR_ALG1).
         */
        template <typename T_Scalar>
        void preprocessSpmm(op_t                             opA,
                            op_t                             opB,
                            const T_Scalar*                  alpha,
                            const cusparse::SpMatBase&       matA,
                            const cusparse::DnMat<T_Scalar>& matB,
                            const T_Scalar*                  beta,
                            cusparse::DnMat<T_Scalar>&       matC,
                            cusparseSpMMAlg_t                alg = cusparseSpMMAlg_t::CUSPARSE_SPMM_CSR_ALG1);

        /**
         * @brief Perform sparse matrix-matrix multiplication.
         * @details Computes C = alpha * op(A) * op(B) + beta * C, where A is a sparse matrix,
         *          B and C are dense matrices, and op(X) is either X, X^T, or X^H.
         * @tparam T_Scalar The scalar type (e.g., float, double, complex).
         * @param[in] opA Operation to apply to matrix A (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] opB Operation to apply to matrix B (e.g., no-transpose, transpose, conjugate-transpose).
         * @param[in] alpha Pointer to scalar alpha.
         * @param[in] matA Sparse matrix A.
         * @param[in] matB Dense matrix B.
         * @param[in] beta Pointer to scalar beta.
         * @param[in,out] matC Output dense matrix C.
         * @param[in] alg Algorithm to use for SpMM (default: CUSPARSE_SPMM_CSR_ALG1).
         */
        template <typename T_Scalar>
        void performSpmm(op_t                             opA,
                         op_t                             opB,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnMat<T_Scalar>& matB,
                         const T_Scalar*                  beta,
                         cusparse::DnMat<T_Scalar>&       matC,
                         cusparseSpMMAlg_t                alg = cusparseSpMMAlg_t::CUSPARSE_SPMM_CSR_ALG1);

    private:
        DeviceBufferVoid& deviceWork() { return m_deviceBuffer; }

    private:
        cusparseHandle_t m_handle{nullptr};
        std::size_t m_workspaceInBytes;
        DeviceBufferVoid m_deviceBuffer;

        void defaults();
};

/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_cusparse_utils
 * @brief Returns the global cuSPARSE handler instance.
 * @details Provides access to a singleton cuSPARSE handler that can be used throughout
 *          the application for cuSPARSE operations.
 * @return Reference to the global @ref CuSparseHandler instance.
 */
inline CuSparseHandler& globalCuSparseHandler()
{
    static CuSparseHandler handler;
    return handler;
}

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUSPARSE_HANDLER_HPP_
