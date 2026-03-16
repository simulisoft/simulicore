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

#ifndef CULITE_LAPACK_BASE_HPP_
#define CULITE_LAPACK_BASE_HPP_

/**
 * @file
 */

#include <culite/types/integer.hpp>
#include <culite/types/enums.hpp>
#include <culite/types/property.hpp>
#include <culite/types/property.hpp>
#include <culite/generic/meta2d.hpp>
#include <culite/generic/cusolver_handler.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

namespace dns { template <typename T_Scalar> class XxVector; }

/**
 * @nosubgrouping
 * @brief The abstract linear solver base for dense device matrices.
 * @details Provides a base class for GPU-based linear solvers that use cuSOLVER decomposition methods.
 *          This is the device version that operates on matrices stored in GPU memory.
 * @tparam T_Matrix The device matrix type (e.g., dns::XxMatrix<T_Scalar>).
 */
template <typename T_Matrix>
class LapackBase {

    private:
        using T_Scalar = typename T_Matrix::value_type;
        using T_Vector = dns::XxVector<T_Scalar>;

    protected:
        LapackBase(CuSolverHandler& cusolver, decomp_t decompType);
        ~LapackBase();

    public:

        /**
         * @brief Reserves workspace memory for matrix decomposition.
         * @details Allocates the necessary device memory buffers required to perform
         *          the decomposition of a matrix with the same dimensions as @p mat.
         * @param[in] mat The matrix whose dimensions determine the workspace size.
         */
        void reserve(const T_Matrix& mat);
        
        /**
         * @brief Clears the solver internal data.
         * @details Clears the solver internal data and resets all settings.
         */
        void clear();
        
        /**
         * @brief Performs matrix decomposition.
         * @details Computes the matrix factorization without modifying the input matrix.
         * @param[in] mat The matrix to be decomposed.
         */
        void decompose(const T_Matrix& mat);
        
        /**
         * @brief Performs in-place matrix solution.
         * @details Solves the linear system using the precomputed factorization.
         * @param[in,out] rhs On input, the right hand side matrix, on exit is overwritten with the solution.
         */
        void solve(T_Matrix& rhs) const;
        
        /**
         * @brief Performs in-place vector solution.
         * @details Solves the linear system using the precomputed factorization.
         * @param[in,out] rhs On input, the right hand side vector, on exit is overwritten with the solution.
         */
        void solve(T_Vector& rhs) const;

    private:

    private:
        CuSolverHandler& m_cusolver;
        const decomp_t m_decompType;
        Meta2D<int_t> m_factorMeta2D;

        decomp_t decompType() const;
        void prepareForDecomposition(const T_Matrix& mat);
        void prepareForSolution(T_Matrix& rhs) const;
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_LAPACK_BASE_HPP_
