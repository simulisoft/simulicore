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
#include "cla3p/lra/lra_math.hpp"

// system

// 3rd

// cla3p
#include "cla3p/lra/lra_xxmatrix.hpp"
#include "cla3p/lra/rank_moderator.hpp"
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/checks/lapack_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/proxies/lapack_proxy.hpp"
#include "cla3p/qr/default_qr.hpp"

#include "cla3p/dense.hpp"
#include "cla3p/algebra.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace lra {
/*-------------------------------------------------*/
template <typename T_Matrix>
static void create_merged_matrices(
		typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
		T_Matrix& mergedLeft, T_Matrix& mergedRight)
{
	int_t m = A.nrows();
	int_t n = A.ncols();
	int_t k1 = A.nrank();
	int_t k2 = B.nrank();
	int_t k = k1 + k2;

	mergedLeft = T_Matrix(m, k);
	mergedRight = T_Matrix(n, k);

	mergedLeft.rblock(0, 0,m,k1) = A.A();
	mergedLeft.rblock(0,k1,m,k2) = B.A();

	mergedLeft.rblock(0, 0,m,k1).iscale(alpha);
	mergedLeft.rblock(0,k1,m,k2).iscale(beta);

	mergedRight.rblock(0, 0,n,k1) = A.B();
	mergedRight.rblock(0,k1,n,k2) = B.B();
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static void update_matrix_with_R_transposed(const DefaultQR<T_Matrix>& qr, T_Matrix& A)
{
	using T_Scalar = typename T_Matrix::value_type;

	const T_Matrix& R = qr.R();

	blas::trmm('R', 'U', 'C', 'N', 
			A.nrows(), A.ncols(), T_Scalar(1), 
			R.values(), R.ld(), 
			A.values(), A.ld());
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static void update_matrix_with_Q(const DefaultQR<T_Matrix>& qr, T_Matrix& A)
{
	using T_Scalar = typename T_Matrix::value_type;
	using T_Vector = dns::XxVector<T_Scalar>;

	const T_Matrix& Q = qr.elementaryReflectors();
	const T_Vector& T = qr.tau();

	int_t info = lapack::xxmqr('L', 'N', 
			A.nrows(), A.ncols(), T.size(),
			Q.values(), A.ld(), 
			T.values(), 
			A.values(), A.ld());

	lapack_info_check(info);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static XxMatrix<T_Matrix> calculate_non_trivial_sum_qq(
		typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
		const RankModerator<T_Matrix>& rmod)
{
	int_t m = A.nrows();
	int_t n = A.ncols();
	//int_t k = A.nrank() + B.nrank(); // k < n && k < m

	T_Matrix mergedLeft;
	T_Matrix mergedRight;
	create_merged_matrices(alpha, A, beta, B, mergedLeft, mergedRight);

	DefaultQR<T_Matrix> qrLeft;
	DefaultQR<T_Matrix> qrRight;
	qrLeft.decompose(mergedLeft);
	qrRight.decompose(mergedRight);

	T_Matrix C = qrLeft.R() * qrRight.R().ctranspose();

	XxMatrix<T_Matrix> tmp = rmod.reduce(C); // (k x k) -> (k x l) &  (k x l)
	XxMatrix<T_Matrix> ret;

	if(tmp) {
		int_t l = tmp.nrank();
		ret = XxMatrix<T_Matrix>(m,n,l);
		ret = 0;
		ret.A().setBlock(0,0,tmp.A());
		ret.B().setBlock(0,0,tmp.B());
		update_matrix_with_Q(qrLeft, ret.A());
		update_matrix_with_Q(qrRight, ret.B());
	} // tmp

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static XxMatrix<T_Matrix> calculate_non_trivial_sum_qb(
		typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
		const RankModerator<T_Matrix>& rmod)
{
	int_t m = A.nrows();
	int_t n = A.ncols();
	//int_t k = A.nrank() + B.nrank(); // n < k < m

	T_Matrix mergedLeft;
	T_Matrix mergedRight;
	create_merged_matrices(alpha, A, beta, B, mergedLeft, mergedRight);

	DefaultQR<T_Matrix> qrLeft;
	qrLeft.decompose(mergedLeft);
	update_matrix_with_R_transposed(qrLeft, mergedRight);

	XxMatrix<T_Matrix> tmp = rmod.reduce(mergedRight); // (n x k) -> (n x l) &  (k x l)
	XxMatrix<T_Matrix> ret;

	if(tmp) {
		int_t l = tmp.nrank();
		ret = XxMatrix<T_Matrix>(m,n,l);
		ret = 0;
		ret.A().setBlock(0,0,tmp.B());
		ret.B().setBlock(0,0,tmp.A());
		update_matrix_with_Q(qrLeft, ret.A());
	} // tmp

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static XxMatrix<T_Matrix> calculate_non_trivial_sum_aq(
		typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
		const RankModerator<T_Matrix>& rmod)
{
	int_t m = A.nrows();
	int_t n = A.ncols();
	//int_t k = A.nrank() + B.nrank(); // m < k < n

	T_Matrix mergedLeft;
	T_Matrix mergedRight;
	create_merged_matrices(alpha, A, beta, B, mergedLeft, mergedRight);

	DefaultQR<T_Matrix> qrRight;
	qrRight.decompose(mergedRight);
	update_matrix_with_R_transposed(qrRight, mergedLeft);

	XxMatrix<T_Matrix> tmp = rmod.reduce(mergedLeft); // (m x k) -> (m x l) &  (k x l)
	XxMatrix<T_Matrix> ret;

	if(tmp) {
		int_t l = tmp.nrank();
		ret = XxMatrix<T_Matrix>(m,n,l);
		ret = 0;
		ret.A().setBlock(0,0,tmp.A());
		ret.B().setBlock(0,0,tmp.B());
		update_matrix_with_Q(qrRight, ret.B());
	} // tmp

	return ret;
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static XxMatrix<T_Matrix> calculate_non_trivial_sum_ab(
		typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
		const RankModerator<T_Matrix>& rmod)
{
	T_Matrix C = alpha * A.A() * A.B().ctranspose() + beta * B.A() * B.B().ctranspose();
	return rmod.reduce(C);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
static XxMatrix<T_Matrix> calculate_non_trivial_sum(
		typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A, 
		typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
		const RankModerator<T_Matrix>& rmod)
{
	int_t m = A.nrows();
	int_t n = A.ncols();
	int_t k = A.nrank() + B.nrank();

	bool qrA = (k < m);
	bool qrB = (k < n);

	/**/ if( qrA &&  qrB) return calculate_non_trivial_sum_qq(alpha, A, beta, B, rmod);
	else if( qrA && !qrB) return calculate_non_trivial_sum_qb(alpha, A, beta, B, rmod);
	else if(!qrA &&  qrB) return calculate_non_trivial_sum_aq(alpha, A, beta, B, rmod);
	else                  return calculate_non_trivial_sum_ab(alpha, A, beta, B, rmod);
}
/*-------------------------------------------------*/
template <typename T_Matrix>
XxMatrix<T_Matrix> add(
    typename T_Matrix::value_type alpha, const XxMatrix<T_Matrix>& A,
    typename T_Matrix::value_type beta, const XxMatrix<T_Matrix>& B,
    const RankModerator<T_Matrix>& rmod)
{
	using T_Scalar = typename T_Matrix::value_type;

	similarity_dim_check(A.nrows(), B.nrows());
	similarity_dim_check(A.ncols(), B.ncols());

	bool validA = alpha != T_Scalar(0) && A.nrank();
	bool validB = beta != T_Scalar(0) && B.nrank();

	XxMatrix<T_Matrix> ret;

	if(validA && validB) {

		ret = calculate_non_trivial_sum(alpha, A, beta, B, rmod);

	} else if(validA && !validB) {

		ret = A;
		ret.iscale(alpha);

	} else if(!validA && validB) {

		ret = B;
		ret.iscale(beta);

	} // cases

	return ret;
}
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
#define instantiate_add(T_Mat) \
template XxMatrix<T_Mat> add( \
    typename T_Mat::value_type, const XxMatrix<T_Mat>&, \
    typename T_Mat::value_type, const XxMatrix<T_Mat>&, \
    const RankModerator<T_Mat>&)
instantiate_add(dns::RdMatrix);
instantiate_add(dns::RfMatrix);
instantiate_add(dns::CdMatrix);
instantiate_add(dns::CfMatrix);
#undef instantiate_add
/*-------------------------------------------------*/
} // namespace lra
} // namespace cla3p
/*-------------------------------------------------*/

