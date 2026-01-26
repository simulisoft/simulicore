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

#ifndef CLA3P_CSX_CHECKS_HPP_
#define CLA3P_CSX_CHECKS_HPP_

/**
 * @file
 */

#include <string>

#include "cla3p/types/property.hpp"
#include "cla3p/checks/basic_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

template <typename T_Int>
void csx_consistency_check(T_Int m, T_Int n, T_Int nnz, 
                           const void *xxxptr, 
                           const void *xxxidx, 
                           const void *values)
{
    if(!(m > 0) || !(n > 0) || !(nnz >= 0)) {
		throw err::NoConsistency("Zero or negative matrix dimensions: (" +
                                 std::to_string(m) + "x" + std::to_string(n) + 
                                 ") with nnz=" + std::to_string(nnz));
	}

	if(!xxxptr) {
		throw err::NoConsistency("Index pointer array is a null pointer.");
	}

	if(nnz && !xxxidx) {
		throw err::NoConsistency("Index array is a null pointer.");
	}

	if(nnz && !values) {
		throw err::NoConsistency("Values array is a null pointer.");
	}
}

/*-------------------------------------------------*/

template <typename T_Int>
void csx_consistency_check(T_Int m, T_Int n, T_Int nnz, 
                           const void *xxxptr, 
                           const void *xxxidx, 
                           const void *values,
                           const Property& prop)
{
	csx_consistency_check(m, n, nnz, xxxptr, xxxidx, values);
	property_compatibility_check(prop, m, n);
}

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_CSX_CHECKS_HPP_
