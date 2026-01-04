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

#ifndef CLA3P_OPERATORS_SCALING_HPP_
#define CLA3P_OPERATORS_SCALING_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

#include "cla3p/virtuals/virtual_object.hpp"
#include "cla3p/virtuals/virtual_scale.hpp"

/*-------------------------------------------------*/

namespace cla3p { 
namespace alias { 

template <typename T_Scalar>
using VirtualScale_vec = VirtualScale<
	dns::XxVector<T_Scalar>,
	VirtualObject<dns::XxVector<T_Scalar>>>;

template <typename T_Scalar>
using VirtualScale_dns = VirtualScale<
	dns::XxMatrix<T_Scalar>,
	VirtualObject<dns::XxMatrix<T_Scalar>>>;

template <typename T_Int, typename T_Scalar>
using VirtualScale_csc = VirtualScale<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>>;

} // namespace alias
} // namespace cla3p
/*-------------------------------------------------*/

/*
 * Generic scale operator
 */
template <typename T_Virtual> 
cla3p::VirtualScale<typename T_Virtual::result_type, T_Virtual>
operator*(
	typename T_Virtual::result_type::value_type val, 
	const cla3p::VirtualExpression<typename T_Virtual::result_type, T_Virtual>& src)
{ 
	return cla3p::VirtualScale<typename T_Virtual::result_type, T_Virtual>(src.self(), val); 
}

/**
 * @ingroup cla3p_module_index_math_operators_scal
 * @brief Multiplies a scalar with a vector.
 * @details Performs the operation <b>val * X</b>
 * @param[in] val The coefficient value.
 * @param[in] X The input vector.
 * @return The virtually scaled vector.
 */
template <typename T_Scalar>
cla3p::alias::VirtualScale_vec<T_Scalar>
operator*(T_Scalar val, const cla3p::dns::XxVector<T_Scalar>& X) 
{ 
	return cla3p::alias::VirtualScale_vec<T_Scalar>(X.virtualize(), val);
}

/**
 * @ingroup cla3p_module_index_math_operators_scal
 * @brief Multiplies a scalar with a dense matrix.
 * @details Performs the operation <b>val * A</b>
 * @param[in] val The coefficient value.
 * @param[in] A The input matrix.
 * @return The virtually scaled matrix.
 */
template <typename T_Scalar>
cla3p::alias::VirtualScale_dns<T_Scalar>
operator*(T_Scalar val, const cla3p::dns::XxMatrix<T_Scalar>& A) 
{ 
	return cla3p::alias::VirtualScale_dns<T_Scalar>(A.virtualize(), val);
}

/**
 * @ingroup cla3p_module_index_math_operators_scal
 * @brief Multiplies a scalar with a sparse matrix.
 * @details Performs the operation <b>val * A</b>
 * @param[in] val The coefficient value.
 * @param[in] A The input matrix.
 * @return The scaled matrix copy.
 */
template <typename T_Int, typename T_Scalar>
cla3p::alias::VirtualScale_csc<T_Int,T_Scalar>
operator*(T_Scalar val, const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A) 
{ 
	return cla3p::alias::VirtualScale_csc<T_Int,T_Scalar>(A.virtualize(), val);
}

/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/

/**
 * @ingroup cla3p_module_index_math_operators_scal
 * @brief Scaling operator.
 * @details Scales `X` by `val`.
 * @param[in] X The input dense vector.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Scalar>
void operator*=(
		cla3p::dns::XxVector<T_Scalar>& X, 
		T_Scalar val)
{
	X.iscale(val);
}

/**
 * @ingroup cla3p_module_index_math_operators_scal
 * @brief Scaling operator.
 * @details Scales `A` by `val`.
 * @param[in] A The input dense matrix.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Scalar>
void operator*=(
		cla3p::dns::XxMatrix<T_Scalar>& A, 
		T_Scalar val)
{
	A.iscale(val);
}

/**
 * @ingroup cla3p_module_index_math_operators_scal
 * @brief Scaling operator.
 * @details Scales `A` by `val`.
 * @param[in] A The input sparse matrix.
 * @param[in] val The scaling coefficient.
 */
template <typename T_Int, typename T_Scalar>
void operator*=(
		cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
		T_Scalar val)
{
	A.iscale(val);
}

/*-------------------------------------------------*/

#endif // CLA3P_OPERATORS_SCALING_HPP_
