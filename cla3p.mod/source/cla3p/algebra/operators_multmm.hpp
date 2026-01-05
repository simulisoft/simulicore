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

#ifndef CLA3P_OPERATORS_MULTMM_HPP_
#define CLA3P_OPERATORS_MULTMM_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/virtuals/virtual_product.hpp"

/*-------------------------------------------------*/

namespace cla3p { 
namespace alias { 

template <typename T_Scalar>
using VirtualProduct_dnsmm = VirtualProduct<
	dns::XxMatrix<T_Scalar>,
	VirtualObject<dns::XxMatrix<T_Scalar>>,
	VirtualObject<dns::XxMatrix<T_Scalar>>>;

template <typename T_Int, typename T_Scalar>
using VirtualProduct_cscmm = VirtualProduct<
	dns::XxMatrix<T_Scalar>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>,
	VirtualObject<dns::XxMatrix<T_Scalar>>>;

template <typename T_Int, typename T_Scalar>
using VirtualProduct_spmm = VirtualProduct<
	csc::XxMatrix<T_Int,T_Scalar>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>,
	VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>>;

} // namespace alias
} // namespace cla3p
/*-------------------------------------------------*/

/*
 * Virtual x XxMatrix
 */

template <typename T_Scalar, typename T_Left>
cla3p::VirtualProduct<
	cla3p::dns::XxMatrix<T_Scalar>,
	T_Left,
	cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>> 
operator*(
	const cla3p::VirtualExpression<typename T_Left::result_type, T_Left>& left, 
	const cla3p::dns::XxMatrix<T_Scalar>& B) 
{ 
	return cla3p::VirtualProduct<
		cla3p::dns::XxMatrix<T_Scalar>,
		T_Left,
		cla3p::VirtualObject<cla3p::dns::XxMatrix<T_Scalar>>>(left.self(), B.virtualize());
}

template <typename T_Int, typename T_Scalar, typename T_Left>
cla3p::VirtualProduct<
	cla3p::csc::XxMatrix<T_Int,T_Scalar>,
	T_Left, 
	cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>
operator*(
	const cla3p::VirtualExpression<typename T_Left::result_type, T_Left>& left, 
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B) 
{ 
	return cla3p::VirtualProduct<
		cla3p::csc::XxMatrix<T_Int,T_Scalar>,
		T_Left,
		cla3p::VirtualObject<cla3p::csc::XxMatrix<T_Int,T_Scalar>>>(left.self(), B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Matrix-matrix multiplication for dense matrices.
 * @details Performs the operation @f$ A \cdot B @f$.
 * @param[in] A The left-hand side input matrix.
 * @param[in] B The right-hand side input matrix.
 * @return The virtual product representing the multiplication.
 */
template <typename T_Scalar>
cla3p::alias::VirtualProduct_dnsmm<T_Scalar>
operator*(
	const cla3p::dns::XxMatrix<T_Scalar>& A, 
	const cla3p::dns::XxMatrix<T_Scalar>& B) 
{
	return cla3p::alias::VirtualProduct_dnsmm<T_Scalar>(A.virtualize(), B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Matrix-matrix multiplication for sparse-dense matrices.
 * @details Performs the operation @f$ A \cdot B @f$ where @p A is sparse and
 *          @p B is dense.
 * @note When forming a virtual product with a sparse matrix and a transposed
 *       dense matrix (e.g., @f$ A \cdot B^T @f$), the library explicitly
 *       calculates the transpose. If used within a loop, pre-calculate the
 *       transpose outside the loop to minimize computational overhead.
 * @param[in] A The left-hand side sparse matrix.
 * @param[in] B The right-hand side dense matrix.
 * @return The resulting dense matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::alias::VirtualProduct_cscmm<T_Int,T_Scalar>
operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::dns::XxMatrix<T_Scalar>& B) 
{
	return cla3p::alias::VirtualProduct_cscmm<T_Int,T_Scalar>(A.virtualize(), B.virtualize());
}

/**
 * @ingroup cla3p_module_index_math_operators_mult
 * @brief Matrix-matrix multiplication for sparse matrices.
 * @details Performs the operation @f$ A \cdot B @f$ where both matrices are sparse.
 * @param[in] A The left-hand side sparse matrix.
 * @param[in] B The right-hand side sparse matrix.
 * @return The resulting sparse matrix.
 */
template <typename T_Int, typename T_Scalar>
cla3p::alias::VirtualProduct_spmm<T_Int,T_Scalar>
operator*(
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& A, 
	const cla3p::csc::XxMatrix<T_Int,T_Scalar>& B) 
{
	return cla3p::alias::VirtualProduct_spmm<T_Int,T_Scalar>(A.virtualize(), B.virtualize());
}

#endif // CLA3P_OPERATORS_MULTMM_HPP_
