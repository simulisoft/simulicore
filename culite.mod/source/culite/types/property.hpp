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

#ifndef CULITE_PROPERTY_HPP_
#define CULITE_PROPERTY_HPP_

/**
 * @file
 */

#include <cla3p/types/property.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_forwards
 * @brief The property class.
 * @details Alias of @ref cla3p::Property into the `culite` namespace.
 *          Used to define matrix properties.
 *          Constructing a matrix with property can accelerate certain operations,
 *          but at the same time imposes equivalent constraints on others.
 */
using Property = ::cla3p::Property;

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_PROPERTY_HPP_
