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

#ifndef CLA3P_LRA_HPP_
#define CLA3P_LRA_HPP_

#include "cla3p/lra/lra_xxmatrix.hpp"
#include "cla3p/lra/rank_moderator.hpp"

namespace cla3p {
namespace lra {

/**
 * @ingroup cla3p_module_index_matrices_lowrank
 * @brief Double precision real low-rank matrix.
 */
using RdMatrix = XxMatrix<dns::RdMatrix>;

/**
 * @ingroup cla3p_module_index_matrices_lowrank
 * @brief Single precision real low-rank matrix.
 */
using RfMatrix = XxMatrix<dns::RfMatrix>;

/**
 * @ingroup cla3p_module_index_matrices_lowrank
 * @brief Double precision complex low-rank matrix.
 */
using CdMatrix = XxMatrix<dns::CdMatrix>;

/**
 * @ingroup cla3p_module_index_matrices_lowrank
 * @brief Single precision complex low-rank matrix.
 */
using CfMatrix = XxMatrix<dns::CfMatrix>;

} // namespace lra
} // namespace cla3p

#endif // CLA3P_LRA_HPP_
