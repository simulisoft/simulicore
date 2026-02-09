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

#ifndef CULITE_VIRTUAL_EXPRESSION_HPP_
#define CULITE_VIRTUAL_EXPRESSION_HPP_

/**
 * @file
 */

/*-------------------------------------------------*/

namespace culite {
namespace dns { template <typename T_Scalar> class XxVector; }
namespace dns { template <typename T_Scalar> class XxMatrix; }
namespace csr { template <typename T_Int, typename T_Scalar> class XxMatrix; }
namespace csc { template <typename T_Int, typename T_Scalar> class XxMatrix; }
} // namespace culite

/*-------------------------------------------------*/
namespace culite { 
/*-------------------------------------------------*/

namespace alias { 

template <typename T_Scalar, typename T_Virtual>
using VirtualExpr_vec = ::cla3p::VirtualExpression<dns::XxVector<T_Scalar>, T_Virtual>;

template <typename T_Scalar, typename T_Virtual>
using VirtualExpr_dns = ::cla3p::VirtualExpression<dns::XxMatrix<T_Scalar>, T_Virtual>;

template <typename T_Int, typename T_Scalar, typename T_Virtual>
using VirtualExpr_csr = ::cla3p::VirtualExpression<csr::XxMatrix<T_Int,T_Scalar>, T_Virtual>;

template <typename T_Int, typename T_Scalar, typename T_Virtual>
using VirtualExpr_csc = ::cla3p::VirtualExpression<csc::XxMatrix<T_Int,T_Scalar>, T_Virtual>;

} // namespace alias

/*-------------------------------------------------*/
} // namespace culite
/*-------------------------------------------------*/

#endif // CULITE_VIRTUAL_EXPRESSION_HPP_
