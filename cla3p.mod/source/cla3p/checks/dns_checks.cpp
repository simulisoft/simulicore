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
#include "cla3p/checks/dns_checks.hpp"

// system

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/error/literals.hpp"
#include "cla3p/checks/basic_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
void dns_consistency_check(int_t m, int_t n, const void *a, int_t lda)
{
	if(!(m > 0) || !(n > 0)) {
		throw err::NoConsistency(msg::InvalidDimensions());
	}

	if(!a) {
		throw err::NoConsistency(msg::InvalidPointer());
	}

	if(lda < m) {
		throw err::NoConsistency(msg::InvalidLeadingDimension());
	}
}
/*-------------------------------------------------*/
void dns_consistency_check(const Property& prop, int_t m, int_t n, const void *a, int_t lda)
{
	dns_consistency_check(m, n, a, lda);
	property_compatibility_check(prop, m, n);
}
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
