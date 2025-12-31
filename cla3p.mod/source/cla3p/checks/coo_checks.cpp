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
#include "cla3p/checks/coo_checks.hpp"

// system

// 3rd

// cla3p
#include "cla3p/support/utils.hpp"
#include "cla3p/error/exceptions.hpp"
#include "cla3p/error/literals.hpp"
#include "cla3p/checks/basic_checks.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/
void coo_consistency_check(const Property& prop, int_t m, int_t n)
{
	if(!m || !n) {
		throw err::NoConsistency(msg::InvalidDimensions());
	}

	property_compatibility_check(prop, m, n);
}
/*-------------------------------------------------*/
void coo_check_coord(int_t m, int_t n, const Property& prop, const int_t& i, const int_t& j)
{
	if(i < 0 || i >= m || j < 0 || j >= n)
	{
		throw err::OutOfBounds(msg::IndexOutOfBounds(m,n,i,j));
	}

	if(prop.isLower() && i < j) {
		throw err::InvalidOp("Invalid insertion of " + coordToString(i,j) + " in upper part");
	}

	if(prop.isUpper() && i > j) {
		throw err::InvalidOp("Invalid insertion of " + coordToString(i,j) + " in lower part");
	}

	if(prop.isSkew() && i == j) {
		throw err::InvalidOp(msg::SkewInconsistency() + " (value on diagonal)");
	}
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void coo_check_triplet(int_t m, int_t n, const Property& prop, const int_t& i, const int_t& j, const T_Scalar& v)
{
	coo_check_coord(m, n, prop, i, j);

	if(prop.isHermitian() && i == j && arith::getIm(v)) {
		throw err::InvalidOp(msg::HermitianInconsistency() + " (imaginary value on diagonal)");
	}
}
/*-------------------------------------------------*/
template void coo_check_triplet(int_t, int_t, const Property&, const int_t&, const int_t&, const real_t&);
template void coo_check_triplet(int_t, int_t, const Property&, const int_t&, const int_t&, const real4_t&);
template void coo_check_triplet(int_t, int_t, const Property&, const int_t&, const int_t&, const complex_t&);
template void coo_check_triplet(int_t, int_t, const Property&, const int_t&, const int_t&, const complex8_t&);
/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/
