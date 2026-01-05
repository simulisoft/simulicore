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

#include "culite/types/traits.hpp"
#include "culite/error/cuda.hpp"
#include "culite/support/imalloc.hpp"
#include "culite/support/utils.hpp"

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
            cusolverStatus_t cusolverStatus = cusolverDnXgetrf_bufferSize(
                handle(),
                params,
                A.nrows(),
                A.ncols(),
                TypeTraits<T_Scalar>::cuda_type(),
                A.values(),
                A.ld(),
                TypeTraits<T_Scalar>::cuda_type(),
                &m_workspaceInBytesOnDevice,
                &m_workspaceInBytesOnHost);
            err::check_cusolver(cusolverStatus);

            ipiv().reserve(std::min(A.nrows(), A.ncols()));
            info().reserve(1);
            factorWork().reserve(A.nrows() * A.ncols() * sizeof(T_Scalar));
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
            cuSolverInt ldf = A.nrows();
            memCopyD2D<T_Scalar>(A.nrows(), A.ncols(), 
                                 A.values(), A.ld(), 
                                 static_cast<T_Scalar*>(factorWork().data()), ldf);
            cusolverDnParams_t params = nullptr;
            cusolverStatus_t cusolverStatus = cusolverDnXgetrf(
                handle(),
                params,
                A.nrows(),
                A.ncols(),
                TypeTraits<T_Scalar>::cuda_type(),
                factorWork().data(),
                ldf,
                ipiv().data(),
                TypeTraits<T_Scalar>::cuda_type(),
                deviceWork().data(),
                m_workspaceInBytesOnDevice,
                hostWork().data(),
                m_workspaceInBytesOnHost,
                info().data());
            err::check_cusolver(cusolverStatus);

            m_factorCudaType = TypeTraits<T_Scalar>::cuda_type();
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
            cuSolverInt ldf = B.nrows();
            cusolverDnParams_t params = nullptr;
            cusolverStatus_t cusolverStatus = cusolverDnXgetrs(
                handle(),
                params,
                CUBLAS_OP_N,
                B.nrows(), // TODO: check dim matching
                B.ncols(),
                m_factorCudaType,
                factorWork().data(),
                ldf,
                ipiv().data(),
                TypeTraits<T_Scalar>::cuda_type(),
                B.values(),
                B.ld(),
                info().data());
            err::check_cusolver(cusolverStatus);
        }

    private:
        DeviceBuffer<cuSolverInt>& ipiv() { return m_ipiv; }
        DeviceBuffer<infoInt>& info() { return m_info; }
        DeviceBufferVoid& factorWork() { return m_factorBuffer; }
        DeviceBufferVoid& deviceWork() { return m_deviceBuffer; }
        DeviceBufferVoid& hostWork() { return m_hostBuffer; }

    private:
        cusolverDnHandle_t m_handle{nullptr};

        size_t m_workspaceInBytesOnDevice;
        size_t m_workspaceInBytesOnHost;
        cudaDataType m_factorCudaType;

        DeviceBuffer<cuSolverInt> m_ipiv;
        DeviceBuffer<infoInt> m_info;
        DeviceBufferVoid m_factorBuffer;
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
