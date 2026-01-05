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

// this file inc
#include "culite/generic/cusolver_handler.hpp"

// system

// 3rd

// culite
#include "culite/error/cuda.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CuSolverHandler::CuSolverHandler()
{
    cusolverStatus_t cusolverStatus = cusolverDnCreate(&m_handle);
    err::check_cusolver(cusolverStatus);
    defaults();
}
/*-------------------------------------------------*/
CuSolverHandler::~CuSolverHandler()
{
    clear();
    cusolverStatus_t cusolverStatus = cusolverDnDestroy(m_handle);
    err::check_cusolver(cusolverStatus);
}
/*-------------------------------------------------*/
void CuSolverHandler::defaults()
{
    m_workspaceInBytesOnDevice = 0;
    m_workspaceInBytesOnHost = 0;
    m_factorCudaType = cudaDataType::CUDA_R_32I; // set to an unsupported type
    //m_factorCudaType = static_cast<cudaDataType>(-1); // unidentified, is this valid?
}
/*-------------------------------------------------*/
void CuSolverHandler::clear()
{
    ipiv().clear();
    info().clear();
    factorWork().clear();
    deviceWork().clear();
    hostWork().clear();

    defaults();
}
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
