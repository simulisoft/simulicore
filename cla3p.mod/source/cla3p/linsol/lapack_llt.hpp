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

#ifndef CLA3P_LAPACK_LLT_HPP_
#define CLA3P_LAPACK_LLT_HPP_

/**
 * @file
 */

#include "cla3p/linsol/lapack_base.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

/**
 * @nosubgrouping
 * @brief The definite Cholesky (LL') linear solver for dense matrices.
 * @details Implements the Cholesky decomposition @f$ A = L L^T @f$ for positive definite matrices.
 * @tparam T_Matrix The matrix type (e.g., dns::XxMatrix<T_Scalar>).
 */
template <typename T_Matrix>
class LapackLLt : public LapackBase<T_Matrix> {

    public:

        // no copy
        LapackLLt(const LapackLLt&) = delete;
        LapackLLt& operator=(const LapackLLt&) = delete;

        /**
         * @brief The default constructor.
         * @details Constructs an empty solver object.
         */
        LapackLLt() : LapackBase<T_Matrix>(decomp_t::LLT) {}

        /**
         * @brief The dimensional constructor.
         * @details Constructs a solver object with @f$ n^2 @f$ preallocated buffered size.
         * @param[in] n The expected problem dimension.
         */
        LapackLLt(int_t n) : LapackBase<T_Matrix>(decomp_t::LLT, n) {}

        /**
         * @brief Destroys the solver.
         * @details Clears all internal data and destroys the solver.
         */
        ~LapackLLt() = default;
};

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_LAPACK_LLT_HPP_
