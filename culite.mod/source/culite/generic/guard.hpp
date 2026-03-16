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

#ifndef CULITE_GUARD_HPP_
#define CULITE_GUARD_HPP_

/**
 * @file
 */

#include <cla3p/generic/guard.hpp>

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_forwards
 * @brief Immutable object wrapper.
 * @tparam T_obj The object type to guard.
 * @details Alias of @ref cla3p::Guard into the `culite` namespace.
 *          Provides a read-only view of an object, preventing modifications
 *          to its contents.
 */
template <typename T_obj>
using Guard = ::cla3p::Guard<T_obj>;

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_GUARD_HPP_
