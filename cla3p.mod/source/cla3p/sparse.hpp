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

#ifndef CLA3P_SPARSE_HPP_
#define CLA3P_SPARSE_HPP_

#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/sparse/coo_xxmatrix.hpp"

namespace cla3p {
namespace csc {

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Double precision real matrix.
 */
using RdMatrix = XxMatrix<int_t,real_t>;

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Single precision real matrix.
 */
using RfMatrix = XxMatrix<int_t,real4_t>;

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Double precision complex matrix.
 */
using CdMatrix = XxMatrix<int_t,complex_t>;

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Single precision complex matrix.
 */
using CfMatrix = XxMatrix<int_t,complex8_t>;

} // namespace csc
} // namespace cla3p


namespace cla3p {
namespace coo {

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Double precision real matrix.
 */
using RdMatrix = XxMatrix<int_t,real_t>;

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Single precision real matrix.
 */
using RfMatrix = XxMatrix<int_t,real4_t>;

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Double precision complex matrix.
 */
using CdMatrix = XxMatrix<int_t,complex_t>;

/**
 * @ingroup cla3p_module_index_matrices_sparse
 * @brief Single precision complex matrix.
 */
using CfMatrix = XxMatrix<int_t,complex8_t>;

} // namespace coo
} // namespace cla3p

#endif // CLA3P_SPARSE_HPP_
