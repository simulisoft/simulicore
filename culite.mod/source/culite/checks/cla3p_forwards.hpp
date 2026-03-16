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

#ifndef CULITE_CHECKS_CLA3P_FORWARDS_HPP_
#define CULITE_CHECKS_CLA3P_FORWARDS_HPP_

/**
 * @file
 */

#include <culite/types/traits.hpp>
#include <culite/types/enums.hpp>
#include <culite/types/property.hpp>

#include <cla3p/checks/basic_checks.hpp>
#include <cla3p/checks/hermitian_coeff_checks.hpp>
#include <cla3p/checks/block_ops_checks.hpp>
#include <cla3p/checks/dns_checks.hpp>
#include <cla3p/checks/csx_checks.hpp>
#include <cla3p/checks/outer_checks.hpp>
#include <cla3p/checks/matrix_math_checks.hpp>
#include <cla3p/checks/decomp_xx_checks.hpp>
#include <cla3p/checks/solve_checks.hpp>

/*-------------------------------------------------*/
namespace culite {
/*-------------------------------------------------*/

template <typename T_Scalar>
void hermitian_coeff_check2(const Property& pr, const T_Scalar& coeff)
{
    using T_Cla3pScalar = typename TypeTraits<T_Scalar>::cla3p_type;
    T_Cla3pScalar cla3pCoeff = TypeTraits<T_Scalar>::toCla3pType(coeff);
    ::cla3p::hermitian_coeff_check<T_Cla3pScalar>(pr, cla3pCoeff);
}

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_CHECKS_CLA3P_FORWARDS_HPP_
