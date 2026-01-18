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
#include "culite/support/utils.hpp"

// system

// 3rd

// culite
#include "culite/error/cuda.hpp"

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/
void syncDevice()
{
	cudaError_t cudaError = cudaDeviceSynchronize();
	err::check_cuda(cudaError);
}
/*-------------------------------------------------*/
void memCopyX2X(std::size_t count, const void *src, void *dest, cudaMemcpyKind kind)
{
	cudaError_t cudaError = cudaMemcpy(dest, src, count, kind);
	err::check_cuda(cudaError);
}
/*-------------------------------------------------*/
void memCopyX2X(std::size_t elemSize, std::size_t m, std::size_t n, 
	            const void *src, std::size_t lds, 
				void *dest, std::size_t ldd, cudaMemcpyKind kind)
{
	size_t copy_width_bytes = m * elemSize; // Height of column in bytes
	size_t src_pitch_bytes  = lds * elemSize;   // Total column stride
	size_t dest_pitch_bytes = ldd * elemSize;  // Total column stride

	cudaError_t cudaError = cudaMemcpy2D(
	    dest,              // Destination
	    dest_pitch_bytes,  // Destination pitch (stride between columns)
	    src,               // Source
	    src_pitch_bytes,   // Source pitch (stride between columns)
	    copy_width_bytes,  // Length of one column to copy (bytes)
	    n,                 // Number of columns to copy
	    kind); 		       // Type of transfer

		err::check_cuda(cudaError);	
}
/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
