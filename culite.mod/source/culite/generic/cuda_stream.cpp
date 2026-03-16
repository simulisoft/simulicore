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
#include <culite/generic/cuda_stream.hpp>

// system

// 3rd

// culite
#include <culite/error/cuda.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
CudaStream::CudaStream()
{
    cudaError_t cudaError = cudaStreamCreate(&m_stream);
    err::check_cuda(cudaError);
}
/*-------------------------------------------------*/
CudaStream::~CudaStream()
{
    cudaError_t cudaError = cudaStreamDestroy(m_stream);
    err::check_cuda(cudaError);
}
/*-------------------------------------------------*/
void CudaStream::sync()
{
    cudaError_t cudaError = cudaStreamSynchronize(m_stream);
    err::check_cuda(cudaError);
}
/*-------------------------------------------------*/
cudaStream_t CudaStream::stream() { return m_stream; }
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
