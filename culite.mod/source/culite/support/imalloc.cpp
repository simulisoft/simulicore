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
#include "culite/support/imalloc.hpp"

// system

// 3rd
#include <cuda_runtime.h>

// culite
#include "culite/error/cuda.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
void* device_alloc(std::size_t size)
{
	void *ret = nullptr;

    if (!size) return ret;

	cudaError_t cudaError = cudaMalloc(&ret, size); err::check_cuda(cudaError);
	
	return ret;
}
/*-------------------------------------------------*/
void device_free(void *ptr) noexcept
{
    cudaFree(ptr);
}
/*-------------------------------------------------*/
void* pinned_alloc(std::size_t size)
{
	void *ret = nullptr;

	if (!size) return ret;

	cudaError_t cudaError = cudaMallocHost(&ret, size); err::check_cuda(cudaError);

	return ret;
}
/*-------------------------------------------------*/
void pinned_free(void *ptr) noexcept
{
    cudaFreeHost(ptr);
}
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
