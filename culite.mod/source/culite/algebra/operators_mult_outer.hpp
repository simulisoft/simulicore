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

#ifndef CULITE_OPERATORS_MULT_OUTER_HPP_
#define CULITE_OPERATORS_MULT_OUTER_HPP_

/**
 * @file
 */

#include "culite/virtuals/virtual_expression.hpp"
#include "culite/virtuals/virtual_object.hpp"
#include "culite/virtuals/virtual_rowvec.hpp"
#include "culite/virtuals/virtual_outer.hpp"

/*-------------------------------------------------*/
namespace culite { 
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; } 
} // namespace culite
/*-------------------------------------------------*/

/*
 * XxVector x VirtualRowvec
 */
template <typename T_Scalar>
culite::VirtualOuter<
    culite::dns::XxMatrix<T_Scalar>,
    culite::alias::VirtualObj_vec<T_Scalar>>
operator*(
        const culite::dns::XxVector<T_Scalar>& x,
        const culite::VirtualRowvec<T_Scalar>& vy)
{
    return culite::VirtualOuter<
        culite::dns::XxMatrix<T_Scalar>,
        culite::alias::VirtualObj_vec<T_Scalar>>(x.virtualize(), vy);
}

/*
 * VirtualExpression x VirtualRowvec
 */
template <typename T_Scalar, typename T_Virtual>

culite::VirtualOuter<culite::dns::XxMatrix<T_Scalar>,T_Virtual>
operator*(
        const culite::alias::VirtualExpr_vec<T_Scalar, T_Virtual>& vx,
        const culite::VirtualRowvec<T_Scalar>& vy)
{
    return culite::VirtualOuter<culite::dns::XxMatrix<T_Scalar>,T_Virtual>(vx, vy);
}

/*-------------------------------------------------*/

#endif // CULITE_OPERATORS_MULT_OUTER_HPP_
