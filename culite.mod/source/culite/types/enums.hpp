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

#ifndef CULITE_ENUMS_HPP_
#define CULITE_ENUMS_HPP_

/** 
 * @file
 */

#include <iostream>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_datatypes
 * @enum alloc_t
 * @brief Specifies the memory allocation strategy and physical residency.
 * @details This enum is used to categorize pointers by their origin, allowing the
 *          runtime to resolve the correct deallocation function (e.g., @c cudaFree vs
 *          @c cudaFreeHost) and optimize data transfer directions.
 */
enum class alloc_t {
	Device       = 0, /**< Device memory */
	Pinned          , /**< Pinned memory */
	Managed         , /**< Managed memory */
	Unregistered      /**< Unknown origin memory */
};

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

/**
 * @ingroup culite_module_index_stream_operators
 * @brief Writes to os the type of prop.
 */
std::ostream& operator<<(std::ostream& os, const culite::alloc_t& prop);

#endif // CULITE_ENUMS_HPP_
