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

#ifndef CULITE_SPARSE_HPP_
#define CULITE_SPARSE_HPP_

#include "culite/sparse/csr_xxmatrix.hpp"
#include "culite/sparse/csc_xxmatrix.hpp"

namespace culite {
namespace csr {

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Double precision real matrix.
 */
using RdMatrix = XxMatrix<int_t,real_t>;

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Single precision real matrix.
 */
using RfMatrix = XxMatrix<int_t,real4_t>;

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Double precision complex matrix.
 */
using CdMatrix = XxMatrix<int_t,complex_t>;

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Single precision complex matrix.
 */
using CfMatrix = XxMatrix<int_t,complex8_t>;

} // namespace csr
} // namespace culite

namespace culite {
namespace csc {

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Double precision real matrix.
 */
using RdMatrix = XxMatrix<int_t,real_t>;

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Single precision real matrix.
 */
using RfMatrix = XxMatrix<int_t,real4_t>;

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Double precision complex matrix.
 */
using CdMatrix = XxMatrix<int_t,complex_t>;

/**
 * @ingroup culite_module_index_matrices_sparse
 * @brief Single precision complex matrix.
 */
using CfMatrix = XxMatrix<int_t,complex8_t>;

} // namespace csc
} // namespace culite

#endif // CULITE_SPARSE_HPP_
