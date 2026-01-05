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
#include "culite/generic/cublas_handler.hpp"

// system

// 3rd

// culite
#include "culite/error/cuda.hpp"
#include "culite/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CuBlasHandler::CuBlasHandler()
{
    cublasStatus_t cublasStatus = cublasCreate(&m_handle);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
CuBlasHandler::~CuBlasHandler()
{
    cublasStatus_t cublasStatus = cublasDestroy(m_handle);
    err::check_cublas(cublasStatus);
}
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
