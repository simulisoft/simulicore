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
void check_cublas(cublasStatus_t status)
{
	if (status != CUBLAS_STATUS_SUCCESS) {
		std::string error_code = std::to_string(status);
		std::string error_name = cublasGetStatusName(status);
		std::string error_details = cublasGetStatusString(status);
		std::string error_msg = "[" + error_code + "]: " + error_name + " - " + error_details;
		throw err::CudaException(error_msg);
	}
}
/*-------------------------------------------------*/
std::string CUSTOMcusolverDnGetErrorName(cusolverStatus_t status)
{
	switch (status) {
		case CUSOLVER_STATUS_SUCCESS: return "CUSOLVER_STATUS_SUCCESS";
		case CUSOLVER_STATUS_NOT_INITIALIZED: return "CUSOLVER_STATUS_NOT_INITIALIZED";
		case CUSOLVER_STATUS_ALLOC_FAILED: return "CUSOLVER_STATUS_ALLOC_FAILED";
		case CUSOLVER_STATUS_INVALID_VALUE: return "CUSOLVER_STATUS_INVALID_VALUE";
		case CUSOLVER_STATUS_ARCH_MISMATCH: return "CUSOLVER_STATUS_ARCH_MISMATCH";
		case CUSOLVER_STATUS_MAPPING_ERROR: return "CUSOLVER_STATUS_MAPPING_ERROR";
		case CUSOLVER_STATUS_EXECUTION_FAILED: return "CUSOLVER_STATUS_EXECUTION_FAILED";
		case CUSOLVER_STATUS_INTERNAL_ERROR: return "CUSOLVER_STATUS_INTERNAL_ERROR";
		case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED: return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
		case CUSOLVER_STATUS_NOT_SUPPORTED: return "CUSOLVER_STATUS_NOT_SUPPORTED";
		case CUSOLVER_STATUS_ZERO_PIVOT: return "CUSOLVER_STATUS_ZERO_PIVOT";
    	case CUSOLVER_STATUS_INVALID_LICENSE: return "CUSOLVER_STATUS_INVALID_LICENSE";
    	case CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED: return "CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED";
    	case CUSOLVER_STATUS_IRS_PARAMS_INVALID: return "CUSOLVER_STATUS_IRS_PARAMS_INVALID";
    	case CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC: return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC";
    	case CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE: return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE";
    	case CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER: return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER";
    	case CUSOLVER_STATUS_IRS_INTERNAL_ERROR: return "CUSOLVER_STATUS_IRS_INTERNAL_ERROR";
    	case CUSOLVER_STATUS_IRS_NOT_SUPPORTED: return "CUSOLVER_STATUS_IRS_NOT_SUPPORTED";
    	case CUSOLVER_STATUS_IRS_OUT_OF_RANGE: return "CUSOLVER_STATUS_IRS_OUT_OF_RANGE";
    	case CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES: return "CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES";
    	case CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED: return "CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED";
    	case CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED: return "CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED";
    	case CUSOLVER_STATUS_IRS_MATRIX_SINGULAR: return "CUSOLVER_STATUS_IRS_MATRIX_SINGULAR";
    	case CUSOLVER_STATUS_INVALID_WORKSPACE: return "CUSOLVER_STATUS_INVALID_WORKSPACE";
		default: return "UNKNOWN_CUSOLVER_STATUS";
	}
}
/*-------------------------------------------------*/
void check_cusolver(cusolverStatus_t status)
{
	if (status != CUSOLVER_STATUS_SUCCESS) {
		std::string error_code = std::to_string(status);
		std::string error_name = CUSTOMcusolverDnGetErrorName(status);
		//std::string error_details = cusolverDnGetErrorString(status); // not available
		std::string error_msg = "[" + error_code + "]: " + error_name; // + " - " + error_details;
		throw err::CudaException(error_msg);
	}
}
/*-------------------------------------------------*/
} // namespace err
} // namespace culite
/*-------------------------------------------------*/
