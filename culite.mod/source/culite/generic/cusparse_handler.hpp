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

#include "culite/support/imalloc.hpp"
#include "culite/proxies/cusparse_proxy.hpp"

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The cuSparse handler class.
 * @details This class provides a wrapper around the cuSPARSE library for performing
 *          sparse linear algebra operations on GPU devices. It manages the cuSPARSE handle
 *          and internal workspace buffers required for factorization and solve operations.
 */
class CuSparseHandler {

    private:
        using cuSparseInt = int64_t;

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
        cusparseHandle_t handle() { return m_handle; }

        /**
         * @brief Clear the internal workspace buffers.
         * @details Releases all internal workspace memory allocated for cuSPARSE operations.
         */
        void clear();

        template <typename T_Scalar>
        void reserveSpmv(::cla3p::op_t                    opA,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnVec<T_Scalar>& vecX,
                         const T_Scalar*                  beta,
                         cusparse::DnVec<T_Scalar>&       vecY,
                         cusparseSpMVAlg_t                alg = cusparseSpMVAlg_t::CUSPARSE_SPMV_CSR_ALG1)
        {
            cusparseStatus_t status = cusparseSpMV_bufferSize(handle(),
                                                              cusparse::cla3pOp2cusparseOp(opA),
                                                              alpha,
                                                              matA.descr(),
                                                              vecX.descr(),
                                                              beta,
                                                              vecY.descr(),
                                                              TypeTraits<T_Scalar>::cuda_type(),
                                                              alg,
                                                              &m_workspaceInBytes);
            err::check_cusparse(status);

            deviceWork().reserve(m_workspaceInBytes);
        }

        template <typename T_Scalar>
        void preprocessSpmv(::cla3p::op_t                    opA,
                            const T_Scalar*                  alpha,
                            const cusparse::SpMatBase&       matA,
                            const cusparse::DnVec<T_Scalar>& vecX,
                            const T_Scalar*                  beta,
                            cusparse::DnVec<T_Scalar>&       vecY,
                            cusparseSpMVAlg_t                alg = cusparseSpMVAlg_t::CUSPARSE_SPMV_CSR_ALG1)
        {
            cusparseStatus_t status = cusparseSpMV_preprocess(handle(),
                                                              cusparse::cla3pOp2cusparseOp(opA),
                                                              alpha,
                                                              matA.descr(),
                                                              vecX.descr(),
                                                              beta,
                                                              vecY.descr(),
                                                              TypeTraits<T_Scalar>::cuda_type(),
                                                              alg,
                                                              deviceWork().data());
            err::check_cusparse(status);
        }

        template <typename T_Scalar>
        void performSpmv(::cla3p::op_t                    opA,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnVec<T_Scalar>& vecX,
                         const T_Scalar*                  beta,
                         cusparse::DnVec<T_Scalar>&       vecY,
                         cusparseSpMVAlg_t                alg = cusparseSpMVAlg_t::CUSPARSE_SPMV_CSR_ALG1)
        {
            cusparseStatus_t status = cusparseSpMV(handle(),
                                                   cusparse::cla3pOp2cusparseOp(opA),
                                                   alpha,
                                                   matA.descr(),
                                                   vecX.descr(),
                                                   beta,
                                                   vecY.descr(),
                                                   TypeTraits<T_Scalar>::cuda_type(),
                                                   alg,
                                                   deviceWork().data());
            err::check_cusparse(status);
        }

        template <typename T_Scalar>
        void reserveSpmm(::cla3p::op_t                    opA,
                         ::cla3p::op_t                    opB,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnMat<T_Scalar>& matB,
                         const T_Scalar*                  beta,
                         cusparse::DnMat<T_Scalar>&       matC,
                         cusparseSpMMAlg_t                alg = cusparseSpMMAlg_t::CUSPARSE_SPMM_CSR_ALG1)
        {
            cusparseStatus_t status = cusparseSpMM_bufferSize(handle(),
                                                              cusparse::cla3pOp2cusparseOp(opA),
                                                              cusparse::cla3pOp2cusparseOp(opB),
                                                              alpha,
                                                              matA.descr(),
                                                              matB.descr(),
                                                              beta,
                                                              matC.descr(),
                                                              TypeTraits<T_Scalar>::cuda_type(),
                                                              alg,
                                                              &m_workspaceInBytes);
            err::check_cusparse(status);

            deviceWork().reserve(m_workspaceInBytes);
        }

        template <typename T_Scalar>
        void preprocessSpmm(::cla3p::op_t                    opA,
                            ::cla3p::op_t                    opB,
                            const T_Scalar*                  alpha,
                            const cusparse::SpMatBase&       matA,
                            const cusparse::DnMat<T_Scalar>& matB,
                            const T_Scalar*                  beta,
                            cusparse::DnMat<T_Scalar>&       matC,
                            cusparseSpMMAlg_t                alg = cusparseSpMMAlg_t::CUSPARSE_SPMM_CSR_ALG1)
        {
            cusparseStatus_t status = cusparseSpMM_preprocess(handle(),
                                                              cusparse::cla3pOp2cusparseOp(opA),
                                                              cusparse::cla3pOp2cusparseOp(opB),
                                                              alpha,
                                                              matA.descr(),
                                                              matB.descr(),
                                                              beta,
                                                              matC.descr(),
                                                              TypeTraits<T_Scalar>::cuda_type(),
                                                              alg,
                                                              deviceWork().data());
            err::check_cusparse(status);
        }



        template <typename T_Scalar>
        void performSpmm(::cla3p::op_t                    opA,
                         ::cla3p::op_t                    opB,
                         const T_Scalar*                  alpha,
                         const cusparse::SpMatBase&       matA,
                         const cusparse::DnMat<T_Scalar>& matB,
                         const T_Scalar*                  beta,
                         cusparse::DnMat<T_Scalar>&       matC,
                         cusparseSpMMAlg_t                alg = cusparseSpMMAlg_t::CUSPARSE_SPMM_CSR_ALG1)
        {
            cusparseStatus_t status = cusparseSpMM(handle(),
                                                   cusparse::cla3pOp2cusparseOp(opA),
                                                   cusparse::cla3pOp2cusparseOp(opB),
                                                   alpha,
                                                   matA.descr(),
                                                   matB.descr(),
                                                   beta,
                                                   matC.descr(),
                                                   TypeTraits<T_Scalar>::cuda_type(),
                                                   alg,
                                                   deviceWork().data());
            err::check_cusparse(status);
        }

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
