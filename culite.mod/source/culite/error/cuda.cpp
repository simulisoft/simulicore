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
#include "culite/error/cuda.hpp"

// system
#include <string>

// 3rd

// culite
#include "culite/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace err {
/*-------------------------------------------------*/
void check_cuda(cudaError_t cudaError)
{
	if (cudaError != cudaSuccess) {
		std::string error_code = std::to_string(cudaError);
		std::string error_name = cudaGetErrorName(cudaError);
		std::string error_details = cudaGetErrorString(cudaError);
		std::string error_msg = "[" + error_code + "]: " + error_name + " - " + error_details;
		throw err::CudaException(error_msg);
	}
}
/*-------------------------------------------------*/
} // namespace err
} // namespace culite
/*-------------------------------------------------*/
