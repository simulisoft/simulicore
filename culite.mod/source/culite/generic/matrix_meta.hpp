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

#ifndef CULITE_MATRIX_META_HPP_
#define CULITE_MATRIX_META_HPP_

/**
 * @file
 */

#include <cla3p/generic/matrix_meta.hpp>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_forwards
 * @brief Matrix metadata container.
 * @tparam T_Int Integer type for dimension values.
 * @details Alias of @ref cla3p::MatrixMeta into the `culite` namespace.
 *          Extends @ref culite::Meta2D with matrix property information, storing
 *          dimensions and structural properties (e.g., symmetric, triangular).
 */
template <typename T_Int>
using MatrixMeta = ::cla3p::MatrixMeta<T_Int>;

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_MATRIX_META_HPP_
