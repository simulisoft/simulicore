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

#ifndef CLA3P_VIRTUAL_PRODUCT_SPECS_IMPL_HPP_
#define CLA3P_VIRTUAL_PRODUCT_SPECS_IMPL_HPP_

/**
 * @file
 */

#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"

#include "cla3p/algebra/functional_update.hpp"
#include "cla3p/algebra/functional_multmv.hpp"
#include "cla3p/algebra/functional_multmm.hpp"

/*-------------------------------------------------*/
namespace cla3p { 
/*-------------------------------------------------*/

//
// General Cases
//

template <typename T_Left, typename T_Right>
void VirtualProductEvaluateOnNewSpec(
	const VirtualExpression<typename T_Left::result_type, T_Left>& left,
	const VirtualExpression<typename T_Right::result_type, T_Right>& right,
	typename T_Right::result_type& dest)
{
	typename T_Left::result_type tmpLeft = left.evaluate();
	typename T_Right::result_type tmpRight = right.evaluate();
	VirtualProduct<
		typename T_Right::result_type, 
		VirtualObject<typename T_Left::result_type>, 
		VirtualObject<typename T_Right::result_type>> tmpProduct(tmpLeft.virtualize(), tmpRight.virtualize());
	tmpProduct.evaluateOnNew(dest);
}

template <typename T_Left, typename T_Right>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualExpression<typename T_Left::result_type, T_Left>& left,
	const VirtualExpression<typename T_Right::result_type, T_Right>& right,
	typename T_Right::result_type& dest)
{
	typename T_Left::result_type tmpLeft = left.evaluate();
	typename T_Right::result_type tmpRight = right.evaluate();
	VirtualProduct<
		typename T_Right::result_type, 
		VirtualObject<typename T_Left::result_type>, 
		VirtualObject<typename T_Right::result_type>> tmpProduct(tmpLeft.virtualize(), tmpRight.virtualize());
	tmpProduct.evaluateOnExisting(dest);
}

template <typename T_Left, typename T_Right>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualExpression<typename T_Left::result_type, T_Left>& left,
	const VirtualExpression<typename T_Right::result_type, T_Right>& right,
	typename T_Right::result_type& dest,
	typename T_Right::result_type::value_type coeff)
{
	typename T_Left::result_type tmpLeft = left.evaluate();
	typename T_Right::result_type tmpRight = right.evaluate();
	VirtualProduct<
		typename T_Right::result_type, 
		VirtualObject<typename T_Left::result_type>, 
		VirtualObject<typename T_Right::result_type>> tmpProduct(tmpLeft.virtualize(), tmpRight.virtualize());
	tmpProduct.accumulateOnExisting(dest, coeff);
}

//
// Scaled General Cases
//

template <typename T_Left, typename T_Right>
void VirtualProductEvaluateOnNewSpec(
	const VirtualScale<typename T_Left::result_type, T_Left>& left,
	const VirtualExpression<typename T_Right::result_type, T_Right>& right,
	typename T_Right::result_type& dest)
{
	VirtualProduct<
		typename T_Right::result_type,
		T_Left, 
		T_Right> tmpProduct(left.get(), right.self());
	tmpProduct.evaluateOnNew(dest);
	dest.iscale(left.coeff());
}

template <typename T_Left, typename T_Right>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualScale<typename T_Left::result_type, T_Left>& left,
	const VirtualExpression<typename T_Right::result_type, T_Right>& right,
	typename T_Right::result_type& dest)
{
	VirtualProduct<
		typename T_Right::result_type,
		T_Left, 
		T_Right> tmpProduct(left.get(), right.self());
	tmpProduct.evaluateOnExisting(dest);
	dest.iscale(left.coeff());
}

template <typename T_Left, typename T_Right>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualScale<typename T_Left::result_type, T_Left>& left,
	const VirtualExpression<typename T_Right::result_type, T_Right>& right,
	typename T_Right::result_type& dest, 
	typename T_Right::result_type::value_type coeff)
{
	VirtualProduct<
		typename T_Right::result_type,
		T_Left, 
		T_Right> tmpProduct(left.get(), right.self());
	tmpProduct.accumulateOnExisting(dest, coeff * left.coeff());
}

//
// Dense Matrix-Vector
//

//
// N x V
//
template <typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	dest = dns::XxVector<T_Scalar>(left.get().nrows());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), op_t::N, left.get(), right.get(), T_Scalar(0), dest);
}

template <typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, op_t::N, left.get(), right.get(), T_Scalar(1), dest);
}

//
// T/C x V
//
template <typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	dest = dns::XxVector<T_Scalar>(left.get().ncols());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), left.op(), left.get(), right.get(), T_Scalar(0), dest);
}

template <typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, left.op(), left.get(), right.get(), T_Scalar(1), dest);
}

//
// Sparse (csc) Matrix-Vector
//

//
// N x V
//
template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	dest = dns::XxVector<T_Scalar>(left.get().nrows());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), op_t::N, left.get(), right.get(), T_Scalar(0), dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, op_t::N, left.get(), right.get(), T_Scalar(1), dest);
}

//
// T/C x V
//
template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	dest = dns::XxVector<T_Scalar>(left.get().ncols());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), left.op(), left.get(), right.get(), T_Scalar(0), dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxVector<T_Scalar>>& right, 
	dns::XxVector<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, left.op(), left.get(), right.get(), T_Scalar(1), dest);
}

//
// (dense Matrix x dense Matrix)

//
// N x N
//
template <typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().nrows(), right.get().ncols());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), op_t::N, left.get(), op_t::N, right.get(), T_Scalar(0), dest);
}

template <typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, op_t::N, left.get(), op_t::N, right.get(), T_Scalar(1), dest);
}

//
// T/C x N
//
template <typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().ncols(), right.get().ncols());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), left.op(), left.get(), op_t::N, right.get(), T_Scalar(0), dest);
}

template <typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, left.op(), left.get(), op_t::N, right.get(), T_Scalar(1), dest);
}

//
// N x T/C
//
template <typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().nrows(), right.get().nrows());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), op_t::N, left.get(), right.op(), right.get(), T_Scalar(0), dest);
}

template <typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, op_t::N, left.get(), right.op(), right.get(), T_Scalar(1), dest);
}

//
// T/C x T/C
//
template <typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().ncols(), right.get().nrows());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), left.op(), left.get(), right.op(), right.get(), T_Scalar(0), dest);
}

template <typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, left.op(), left.get(), right.op(), right.get(), T_Scalar(1), dest);
}

//
// (csc Matrix x dense Matrix) Remaining Cases in .cpp
//

//
// N x N
//
template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().nrows(), right.get().ncols());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), op_t::N, left.get(), right.get(), T_Scalar(0), dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, op_t::N, left.get(), right.get(), T_Scalar(1), dest);
}

//
// T/C x N
//
template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().ncols(), right.get().ncols());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), left.op(), left.get(), right.get(), T_Scalar(0), dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, left.op(), left.get(), right.get(), T_Scalar(1), dest);
}

//
// N x T/C
//
template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().nrows(), right.get().nrows());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), op_t::N, left.get(), right.evaluate(), T_Scalar(0), dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, op_t::N, left.get(), right.evaluate(), T_Scalar(1), dest);
}

//
// T/C x T/C
//
template <typename T_Int,typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	dest = dns::XxMatrix<T_Scalar>(left.get().ncols(), right.get().nrows());
	VirtualProductEvaluateOnExistingSpec(left, right, dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest)
{
	ops::mult(T_Scalar(1), left.op(), left.get(), right.evaluate(), T_Scalar(0), dest);
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<dns::XxMatrix<T_Scalar>>& right, 
	dns::XxMatrix<T_Scalar>& dest,
	T_Scalar coeff)
{
	ops::mult(coeff, left.op(), left.get(), right.evaluate(), T_Scalar(1), dest);
}

//
// (csc Matrix x csc Matrix)
//

//
// N x N
//

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest.clear();
	dest = ops::mult(T_Scalar(1), op_t::N, left.get(), op_t::N, right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest = ops::mult(T_Scalar(1), op_t::N, left.get(), op_t::N, right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest,
	T_Scalar coeff)
{
	csc::XxMatrix<T_Int,T_Scalar> tmp = ops::mult(T_Scalar(1), op_t::N, left.get(), op_t::N, right.get());
	ops::update(coeff, tmp, dest);
}

//
// T/C x N
//

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest.clear();
	dest = ops::mult(T_Scalar(1), left.op(), left.get(), op_t::N, right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest = ops::mult(T_Scalar(1), left.op(), left.get(), op_t::N, right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest,
	T_Scalar coeff)
{
	csc::XxMatrix<T_Int,T_Scalar> tmp = ops::mult(T_Scalar(1), left.op(), left.get(), op_t::N, right.get());
	ops::update(coeff, tmp, dest);
}

//
// N x T/C
//

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest.clear();
	dest = ops::mult(T_Scalar(1), op_t::N, left.get(), right.op(), right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest = ops::mult(T_Scalar(1), op_t::N, left.get(), right.op(), right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualObject<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest,
	T_Scalar coeff)
{
	csc::XxMatrix<T_Int,T_Scalar> tmp = ops::mult(T_Scalar(1), op_t::N, left.get(), right.op(), right.get());
	ops::update(coeff, tmp, dest);
}

//
// T/C x T/C
//

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnNewSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest.clear();
	dest = ops::mult(T_Scalar(1), left.op(), left.get(), right.op(), right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductEvaluateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest)
{
	dest = ops::mult(T_Scalar(1), left.op(), left.get(), right.op(), right.get());
}

template <typename T_Int, typename T_Scalar>
void VirtualProductAccumulateOnExistingSpec(
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& left, 
	const VirtualTranspose<csc::XxMatrix<T_Int,T_Scalar>>& right, 
	csc::XxMatrix<T_Int,T_Scalar>& dest,
	T_Scalar coeff)
{
	csc::XxMatrix<T_Int,T_Scalar> tmp = ops::mult(T_Scalar(1), left.op(), left.get(), right.op(), right.get());
	ops::update(coeff, tmp, dest);
}

/*-------------------------------------------------*/
} // namespace cla3p
/*-------------------------------------------------*/

#endif // CLA3P_VIRTUAL_PRODUCT_SPECS_IMPL_HPP_
