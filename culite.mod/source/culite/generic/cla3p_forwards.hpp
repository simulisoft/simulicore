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

#ifndef CULITE_GENERIC_CLA3P_FORWARDS_HPP_
#define CULITE_GENERIC_CLA3P_FORWARDS_HPP_

/**
 * @file
 */

#include <cla3p/generic/ownership.hpp>
#include <cla3p/generic/guard.hpp>
#include <cla3p/generic/meta1d.hpp>
#include <cla3p/generic/meta2d.hpp>
#include <cla3p/generic/matrix_meta.hpp>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @brief Immutable object wrapper.
 * @tparam T_obj The object type to guard.
 * @details Alias of @ref cla3p::Guard into the `culite` namespace.
 *          Provides a read-only view of an object, preventing modifications
 *          to its contents.
 */
template <typename T_obj>
using Guard = ::cla3p::Guard<T_obj>;

/**
 * @brief One-dimensional metadata container.
 * @tparam T_Int Integer type for storing the dimension value.
 * @details Alias of @ref cla3p::Meta1D into the `culite` namespace.
 *          Provides storage and accessors for vector size.
 */
template <typename T_Int>
using Meta1D = ::cla3p::Meta1D<T_Int>;

/**
 * @brief Two-dimensional metadata container.
 * @tparam T_Int Integer type for storing dimension values.
 * @details Alias of @ref cla3p::Meta2D into the `culite` namespace.
 *          Provides storage and accessors for matrix dimensions (rows and columns).
 */
template <typename T_Int>
using Meta2D = ::cla3p::Meta2D<T_Int>;

/**
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

#endif // CULITE_GENERIC_CLA3P_FORWARDS_HPP_
