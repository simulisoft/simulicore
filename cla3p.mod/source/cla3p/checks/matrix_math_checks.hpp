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

#ifndef CLA3P_MATRIX_MATH_CHECKS_HPP_
#define CLA3P_MATRIX_MATH_CHECKS_HPP_

/**
 * @file
 */

#include <string>
#include <sstream>

#include "cla3p/types/enums.hpp"
#include "cla3p/types/integer.hpp"
#include "cla3p/types/property.hpp"
#include "cla3p/types/operation.hpp"
#include "cla3p/generic/meta1d.hpp"
#include "cla3p/generic/meta2d.hpp"
#include "cla3p/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

template <typename T_Int>
void mult_dim_check(op_t opA, const Meta2D<T_Int>& metaA,
                    op_t opB, const Meta2D<T_Int>& metaB,
                    const Meta2D<T_Int>& metaC)
{
	T_Int m = (opA == op_t::N ? metaA.nrows() : metaA.ncols());
	T_Int n = (opB == op_t::N ? metaB.ncols() : metaB.nrows());

	T_Int kA = (opA == op_t::N ? metaA.ncols() : metaA.nrows());
	T_Int kB = (opB == op_t::N ? metaB.nrows() : metaB.ncols());

	if(metaC.nrows() != m || metaC.ncols() != n || kA != kB) {
        std::stringstream ss;
        ss << "Dimension mismatch for matrix multiplication: \n"
           << "A: " << metaA << " opA: " << opA << "\n"
           << "B: " << metaB << " opB: " << opB << "\n"
           << "C: " << metaC;
		throw err::NoConsistency(ss.str());
	}
}

template <typename T_Int>
void mult_dim_check(op_t opA, const Meta2D<T_Int>& metaA,
                    const Meta1D<T_Int>& metaX,
                    const Meta1D<T_Int>& metaY)
{
	T_Int m = (opA == op_t::N ? metaA.nrows() : metaA.ncols());
	T_Int n = (opA == op_t::N ? metaA.ncols() : metaA.nrows());

	if(metaX.size() != m || metaY.size() != n) {
        std::stringstream ss;
        ss << "Dimension mismatch for matrix-vector multiplication: \n"
           << "A: " << metaA << " opA: " << opA << "\n"
           << "x: " << metaX << "\n"
           << "y: " << metaY;
        throw err::NoConsistency(ss.str());
	}
}

/*-------------------------------------------------*/

void mult_dim_check(
		int_t nrowsA, int_t ncolsA, const Operation& opA, 
		int_t nrowsB, int_t ncolsB, const Operation& opB, 
		int_t nrowsC, int_t ncolsC);

void mat_x_vec_mult_check(const Operation& opA, 
		const Property& prA, int_t nrowsA, int_t ncolsA, 
		int_t sizeX, int_t sizeY);

#if 0
void mat_x_mat_mult_check(
		const Property& prA, int_t nrowsA, int_t ncolsA, const Operation& opA, 
		const Property& prB, int_t nrowsB, int_t ncolsB, const Operation& opB, 
		const Property& prC, int_t nrowsC, int_t ncolsC);
#endif

void trivec_mult_replace_check(const Property& prA, 
		int_t nrowsA, int_t ncolsA, const Operation& opA, 
		int_t sizeX);

void trimat_mult_replace_check(side_t sideA, 
		const Property& prA, int_t nrowsA, int_t ncolsA, const Operation& opA, 
		const Property& prB, int_t nrowsB, int_t ncolsB);

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_MATRIX_MATH_CHECKS_HPP_
