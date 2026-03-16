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

#ifndef CULITE_META1D_HPP_
#define CULITE_META1D_HPP_

/**
 * @file
 */

#include <cla3p/generic/meta1d.hpp>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_forwards
 * @brief One-dimensional metadata container.
 * @tparam T_Int Integer type for storing the dimension value.
 * @details Alias of @ref cla3p::Meta1D into the `culite` namespace.
 *          Provides storage and accessors for vector size.
 */
template <typename T_Int>
using Meta1D = ::cla3p::Meta1D<T_Int>;

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_META1D_HPP_
