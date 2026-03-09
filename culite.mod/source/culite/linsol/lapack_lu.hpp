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

#ifndef CULITE_LAPACK_LU_HPP_
#define CULITE_LAPACK_LU_HPP_

/**
 * @file
 */

#include "culite/linsol/lapack_base.hpp"

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }

/**
 * @nosubgrouping
 * @brief The partial pivoting LU linear solver for dense device matrices.
 * @details Implements the LU decomposition with partial pivoting @f$ P A = L U @f$ for general matrices
 *          stored in device memory. This is the GPU version using cuSOLVER.
 * @tparam T_Matrix The device matrix type (e.g., dns::XxMatrix<T_Scalar>).
 */
template <typename T_Matrix>
class LapackLU : public LapackBase<T_Matrix> {

    public:

        // no copy
        LapackLU(const LapackLU&) = delete;
        LapackLU& operator=(const LapackLU&) = delete;

    /**
     * @brief Constructor.
     * @details Initializes the LU solver with the specified cuSOLVER handler.
     * @param[in] cusolver Reference to the cuSOLVER handler (defaults to global handler).
     */
    LapackLU(CuSolverHandler& cusolver = globalCuSolverHandler()) 
    : LapackBase<T_Matrix>(cusolver, ::cla3p::decomp_t::LU) {}

    /**
     * @brief Destructor.
     * @details Destroys the solver and releases allocated resources.
     */
        ~LapackLU() = default;
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_LAPACK_LU_HPP_