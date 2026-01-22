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

#include <cusparse.h>

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

    private:
        cusparseHandle_t m_handle{nullptr};
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
