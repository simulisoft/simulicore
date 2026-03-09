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

#ifndef CLA3P_OUTER_CHECKS_HPP_
#define CLA3P_OUTER_CHECKS_HPP_

/**
 * @file
 */

#include <sstream>

#include "cla3p/types/property.hpp"
#include "cla3p/generic/meta1d.hpp"
#include "cla3p/generic/meta2d.hpp"
#include "cla3p/generic/matrix_meta.hpp"
#include "cla3p/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace cla3p {
/*-------------------------------------------------*/

template <typename T_Int>
void outer_product_consistency_check(bool conjop,
                                     const MatrixMeta<T_Int>& metaA, 
                                     const Meta1D<T_Int>& metaX, 
                                     const Meta1D<T_Int>& metaY)
{
    if(!(metaA.prop().isGeneral() || metaA.prop().isSymmetric() || metaA.prop().isHermitian())) {
        std::stringstream ss;
        ss << "Invalid matrix property for outer product: " << metaA.prop();
        throw err::NoConsistency(ss.str());
    } // valid props

    if (metaA.prop().isSymmetric() && conjop) {
        throw err::NoConsistency("For Symmetric matrices, the transpose outer product must be used.");
    }

    if(metaA.prop().isHermitian() && !conjop) {
        throw err::NoConsistency("For Hermitian matrices, the conjugate outer product must be used.");
    }

    if(metaA.nrows() != metaX.size() || metaA.ncols() != metaY.size()) {
        std::stringstream ss;
        ss << "Inconsistent dimensions for outer product: A" << metaA << " x" << metaX << ", y" << metaY;
        throw err::NoConsistency(ss.str());
    }
}

/*-------------------------------------------------*/

template <typename T_Int>
void outer_product_consistency_check(bool conjop, 
                                     T_Int nrowsA, 
                                     T_Int ncolsA, 
                                     const Property& prA, 
                                     T_Int sizeX, 
                                     T_Int sizeY)
{
    MatrixMeta<T_Int> metaA(nrowsA, ncolsA, prA);
    Meta1D<T_Int> metaX(sizeX);
    Meta1D<T_Int> metaY(sizeY);
    outer_product_consistency_check(conjop, metaA, metaX, metaY);
}

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_OUTER_CHECKS_HPP_
