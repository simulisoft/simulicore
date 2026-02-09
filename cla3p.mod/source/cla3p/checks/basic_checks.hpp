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

#ifndef CLA3P_BASIC_CHECKS_HPP_
#define CLA3P_BASIC_CHECKS_HPP_

/**
 * @file
 */
#include <string>

#include "cla3p/types/enums.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

template <typename T_Int>
void square_check(const T_Int& m, const T_Int& n)
{
	if(m != n) {
		std::string errMsg = "Matrix is not square: (" + std::to_string(m) + "x" + std::to_string(n) + ")";
		throw err::NoConsistency(errMsg);
	}
}

template <typename T_Int>
void square_check(const Meta2D<T_Int>& meta)
{
	square_check(meta.nrows(), meta.ncols());
}

/*-------------------------------------------------*/

template <typename T_Int>
void property_compatibility_check(const Property& prop, const T_Int& m, const T_Int& n)
{
	if(!prop.isValid()) {
		std::string errMsg = "Invalid property.";
		throw err::NoConsistency(errMsg);
	}

	if(prop.isSquare()) {
		square_check(m, n);
	}
}

template <typename T_Int>
void property_compatibility_check(const MatrixMeta<T_Int>& meta)
{
	property_compatibility_check(meta.prop(), meta.nrows(), meta.ncols());
}

/*-------------------------------------------------*/

template <typename T_Int>
void similarity_dim_check(const T_Int& n1, const T_Int& n2)
{
	if(n1 != n2) {
		std::string errMsg = "Mismatching dimensions: (" + std::to_string(n1) + " vs " + std::to_string(n2) + ")";
		throw err::NoConsistency(errMsg);
	}
}

/*-------------------------------------------------*/

template <typename T_Int>
void similarity_check(
		const Property& prop1, const T_Int& nrows1, const T_Int& ncols1, 
		const Property& prop2, const T_Int& nrows2, const T_Int& ncols2)
{
	similarity_dim_check(nrows1, nrows2);
	similarity_dim_check(ncols1, ncols2);

	if(prop1 != prop2) {
		std::string errMsg = "Mismatching properties " + prop1.name() + " vs " + prop2.name() + ".";
		throw err::NoConsistency(errMsg);
	}
}

template <typename T_Int>
void similarity_check(const MatrixMeta<T_Int>& meta1,
					  const MatrixMeta<T_Int>& meta2)
{
	similarity_check(meta1.prop(), meta1.nrows(), meta1.ncols(),
					 meta2.prop(), meta2.nrows(), meta2.ncols());
}

template <typename T_Int>
void similarity_check(op_t op1, const MatrixMeta<T_Int>& meta1,
					  op_t op2, const MatrixMeta<T_Int>& meta2)
{
    T_Int m1 = (op1 == op_t::N ? meta1.nrows() : meta1.ncols());
    T_Int n1 = (op1 == op_t::N ? meta1.ncols() : meta1.nrows());

    T_Int m2 = (op2 == op_t::N ? meta2.nrows() : meta2.ncols());
    T_Int n2 = (op2 == op_t::N ? meta2.ncols() : meta2.nrows());

	similarity_check(meta1.prop(), m1, n1,
					 meta2.prop(), m2, n2);
}

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_BASIC_CHECKS_HPP_
