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
#include "cla3p/algebra/functional_add.hpp"

// system

// 3rd

// cla3p
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/bulk/dns_math.hpp"
#include "cla3p/bulk/csc_math.hpp"
#include "cla3p/lra/lra_math.hpp"

#include "cla3p/dense.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
dns::XxVector<T_Scalar> add(
		T_Scalar alpha, const dns::XxVector<T_Scalar>& X,
		T_Scalar beta , const dns::XxVector<T_Scalar>& Y)
{
	similarity_dim_check(X.size(), Y.size());

	dns::XxVector<T_Scalar> ret(X.size());

	blk::dns::add(uplo_t::Full, X.size(), 1, 
		alpha, X.values(), X.size(),
		beta , Y.values(), Y.size(), 
		ret.values(), ret.size());

	return ret;
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template dns::XxVector<T_Scl> add( \
		T_Scl, const dns::XxVector<T_Scl>&, \
		T_Scl, const dns::XxVector<T_Scl>&)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
template <typename T_Scalar>
 dns::XxMatrix<T_Scalar> add(
		T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A,
		T_Scalar beta , const dns::XxMatrix<T_Scalar>& B)
{
	similarity_check(
			A.prop(), A.nrows(), A.ncols(),
			B.prop(), B.nrows(), B.ncols());

	int_t m = A.nrows();
	int_t n = A.ncols();
	Property pr = A.prop();

	dns::XxMatrix<T_Scalar> ret(m, n, pr);

	blk::dns::add(pr.uplo(), m, n, 
			alpha, A.values(), A.ld(),
			beta, B.values(), B.ld(), 
			ret.values(), ret.ld());

	return ret;
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template dns::XxMatrix<T_Scl> add( \
		T_Scl, const dns::XxMatrix<T_Scl>&, \
		T_Scl, const dns::XxMatrix<T_Scl>&)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
csc::XxMatrix<T_Int,T_Scalar> add(
		T_Scalar alpha, const csc::XxMatrix<T_Int,T_Scalar>& A,
		T_Scalar beta, const csc::XxMatrix<T_Int,T_Scalar>& B)
{
	similarity_check(
			A.prop(), A.nrows(), A.ncols(),
			B.prop(), B.nrows(), B.ncols());

	csc::XxMatrix<T_Int,T_Scalar> ret;

	if(beta == T_Scalar(0) && alpha != T_Scalar(0)) {
		ret = A;
		ret.iscale(alpha);
		return ret;
	} else if(alpha == T_Scalar(0) && beta != T_Scalar(0)) {
		ret = B;
		ret.iscale(beta);
		return ret;
	} else if(alpha == T_Scalar(0) && beta == T_Scalar(0)) {
		return ret;
	} // alpha/beta

	int_t    nrowsC  = A.nrows();
	int_t    ncolsC  = A.ncols();
	Property  propC   = A.prop();
	int_t*    colptrC = nullptr;
	int_t*    rowidxC = nullptr;
	T_Scalar* valuesC = nullptr;

	blk::csc::add(nrowsC, ncolsC, 
			alpha, A.colptr(), A.rowidx(), A.values(),
			beta, B.colptr(), B.rowidx(), B.values(),
			&colptrC, &rowidxC, &valuesC);

	ret = csc::XxMatrix<T_Int,T_Scalar>(nrowsC, ncolsC, colptrC, rowidxC, valuesC, true, propC);

	return ret;
}
/*-------------------------------------------------*/
#define instantiate_add(T_Int,T_Scl) \
template csc::XxMatrix<T_Int,T_Scl> add( \
		T_Scl, const csc::XxMatrix<T_Int,T_Scl>&, \
		T_Scl, const csc::XxMatrix<T_Int,T_Scl>&)
instantiate_add(int_t,real_t);
instantiate_add(int_t,real4_t);
instantiate_add(int_t,complex_t);
instantiate_add(int_t,complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
template <typename T_Matrix>
lra::XxMatrix<T_Matrix> add(
		typename T_Matrix::value_type alpha, const lra::XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const lra::XxMatrix<T_Matrix>& B,
		const lra::RankModerator<T_Matrix>& rmod)
{
	return lra::add(alpha, A, beta, B, rmod);
}
/*-------------------------------------------------*/
#define instantiate_add(T_Mat) \
template lra::XxMatrix<T_Mat> add( \
		typename T_Mat::value_type, const lra::XxMatrix<T_Mat>&, \
		typename T_Mat::value_type, const lra::XxMatrix<T_Mat>&, \
		const lra::RankModerator<T_Mat>&)
instantiate_add(dns::RdMatrix);
instantiate_add(dns::RfMatrix);
instantiate_add(dns::CdMatrix);
instantiate_add(dns::CfMatrix);
#undef instantiate_add
/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/
