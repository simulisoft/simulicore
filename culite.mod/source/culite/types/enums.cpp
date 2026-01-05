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

// this file inc
#include "culite/types/enums.hpp"

// system

// 3rd
#include <cla3p/types/literals.hpp>

// culite

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const culite::alloc_t& alloc_type)
{
	std::string str = ::cla3p::msg::Unknown();

	if(alloc_type == culite::alloc_t::Device      ) str = "Device Memory";
	if(alloc_type == culite::alloc_t::Pinned      ) str = "Pinned Memory";
	if(alloc_type == culite::alloc_t::Managed     ) str = "Unified Memory";
	if(alloc_type == culite::alloc_t::Unregistered) str = "Unknown Origin Memory";

	os << str;
	return os;
}
/*-------------------------------------------------*/
