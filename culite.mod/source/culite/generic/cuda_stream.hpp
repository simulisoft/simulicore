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

#ifndef CCULITE_CUDA_STREAM_HPP_
#define CCULITE_CUDA_STREAM_HPP_

/**
 * @file
 */

#include <cuda_runtime.h>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping 
 * @brief CUDA stream encapsulation class.
 * @details This class provides a RAII wrapper around CUDA streams for managing
 *          asynchronous kernel execution and memory operations. The stream is
 *          created on construction and destroyed on destruction.
 */
class CudaStream {

    public:

        // disable copying
        CudaStream(const CudaStream&) = delete;
        CudaStream& operator=(const CudaStream&) = delete;
        
        /**
         * @brief Constructor.
         * @details Creates a new CUDA stream for asynchronous operations.
         */
        CudaStream();

        /**
         * @brief Destructor.
         * @details Destroys the CUDA stream and releases associated resources.
         */
        ~CudaStream();
        
        /**
         * @brief Synchronizes the CUDA stream.
         * @details Blocks the host until all operations in the stream are complete.
         */
        void sync();

        /**
         * @brief Retrieves the underlying CUDA stream handle.
         * @details Provides access to the native CUDA stream handle for use with
         *          CUDA runtime API functions and cuBLAS operations.
         * @return The CUDA stream handle.
         */
        cudaStream_t stream();

    private:
        cudaStream_t m_stream{nullptr};
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CUDA_STREAM_HPP_
