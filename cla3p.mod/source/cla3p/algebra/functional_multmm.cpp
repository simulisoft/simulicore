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
#include "cla3p/algebra/functional_multmm.hpp"

// system
#include <string>

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/checks/matrix_math_checks.hpp"
#include "cla3p/checks/hermitian_coeff_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/bulk/dns_math.hpp"
#include "cla3p/bulk/csc_math.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/algebra/functional_update.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Rhs, typename T_Lhs, typename T_Ret>
static void throw_prop_compatibility_error(
		const T_Rhs& A, 
		const T_Lhs& B, 
		const T_Ret& C)
{
	std::string message = "Bad matrix property combo:\n";
	message = message + "A: " + A.prop().name() + "\n";
	message = message + "B: " + B.prop().name() + "\n";
	message = message + "C: " + C.prop().name();

	throw err::InvalidOp(message);
}
/*-------------------------------------------------*/
template <typename T_Rhs, typename T_Lhs>
static void throw_prop_compatibility_error(
		const T_Rhs& A, 
		const T_Lhs& B) 
{
	std::string message = "Bad matrix property combo:\n";
	message = message + "A: " + A.prop().name() + "\n";
	message = message + "B: " + B.prop().name();

	throw err::InvalidOp(message);
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void mult(T_Scalar alpha,
    op_t opA, const dns::XxMatrix<T_Scalar>& A,
    op_t opB, const dns::XxMatrix<T_Scalar>& B,
		T_Scalar beta, dns::XxMatrix<T_Scalar>& C)
{
	if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = op_t::N;
	if(B.prop().isSymmetric() || B.prop().isHermitian()) opB = op_t::N;

	opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);
	opB = (TypeTraits<T_Scalar>::is_real() && opB == op_t::C ? op_t::T : opB);

	Operation _opA(opA);
	Operation _opB(opB);

	mult_dim_check(
			A.nrows(), A.ncols(), _opA, 
			B.nrows(), B.ncols(), _opB, 
			C.nrows(), C.ncols());

	hermitian_coeff_check(C.prop(), alpha);
	hermitian_coeff_check(C.prop(), beta);

	if(A.prop().isGeneral() && B.prop().isGeneral()) {

		int_t k = (_opA.isTranspose() ? A.nrows() : A.ncols());

		if(C.prop().isGeneral()) {

			blk::dns::gem_x_gem(
					C.nrows(), 
					C.ncols(), 
					k, alpha, 
					opA, A.values(), A.ld(), 
					opB, B.values(), B.ld(), 
					beta, C.values(), C.ld());

		} else if(C.prop().isSymmetric() || C.prop().isHermitian()) {

			blas::gemmt(C.prop().cuplo(), 
					_opA.ctype(), 
					_opB.ctype(), 
					C.ncols(), k, alpha, 
					A.values(), A.ld(), 
					B.values(), B.ld(), 
					beta, 
					C.values(), C.ld());

		} else {

			throw_prop_compatibility_error(A, B, C);

		} // property combos

	} else if(A.prop().isSymmetric() && B.prop().isGeneral() && C.prop().isGeneral()) {

		blk::dns::sym_x_gem(A.prop().uplo(), 
				C.nrows(), 
				C.ncols(), 
				alpha, 
				A.values(), A.ld(), 
				B.values(), B.ld(), 
				beta, 
				C.values(), C.ld());

	} else if(A.prop().isHermitian() && B.prop().isGeneral() && C.prop().isGeneral()) {

		blk::dns::hem_x_gem(
				A.prop().uplo(), 
				C.nrows(), 
				C.ncols(), 
				alpha, 
				A.values(), A.ld(), 
				B.values(), B.ld(), 
				beta, 
				C.values(), C.ld());

	} else if(A.prop().isGeneral() && B.prop().isSymmetric() && C.prop().isGeneral()) {

		blk::dns::gem_x_sym(B.prop().uplo(), 
				C.nrows(), 
				C.ncols(), 
				alpha, 
				B.values(), B.ld(), 
				A.values(), A.ld(), 
				beta, 
				C.values(), C.ld());

	} else if(A.prop().isGeneral() && B.prop().isHermitian() && C.prop().isGeneral()) {

		blk::dns::gem_x_hem(B.prop().uplo(), 
				C.nrows(), 
				C.ncols(), 
				alpha, 
				B.values(), B.ld(), 
				A.values(), A.ld(), 
				beta, 
				C.values(), C.ld());

	} else if(A.prop().isTriangular() && B.prop().isGeneral() && C.prop().isGeneral()) {

		if(beta == T_Scalar(0)) {

			blk::dns::trm_x_gem(A.prop().uplo(), opA, 
					C.nrows(), 
					C.ncols(), 
					B.nrows(), 
					alpha, 
					A.values(), A.ld(), 
					B.values(), B.ld(), 
					C.values(), C.ld());

		} else {

			dns::XxMatrix<T_Scalar> tmp(C.nrows(), C.ncols());
			blk::dns::trm_x_gem(A.prop().uplo(), opA, 
					C.nrows(), 
					C.ncols(), 
					B.nrows(), 
					alpha, 
					A.values(), A.ld(), 
					B.values(), B.ld(), 
					tmp.values(), tmp.ld());
			C.iscale(beta);
			ops::update(T_Scalar(1), tmp, C);

		} // beta

	} else if(A.prop().isGeneral() && B.prop().isTriangular() && C.prop().isGeneral()) {

		if(beta == T_Scalar(0)) {

			blk::dns::gem_x_trm(B.prop().uplo(), opB, 
					C.nrows(), 
					C.ncols(), 
					A.nrows(), 
					alpha, 
					B.values(), B.ld(), 
					A.values(), A.ld(), 
					C.values(), C.ld());

		} else {

			dns::XxMatrix<T_Scalar> tmp(C.nrows(), C.ncols());
			blk::dns::gem_x_trm(B.prop().uplo(), opB, 
					C.nrows(), 
					C.ncols(), 
					A.nrows(), 
					alpha, 
					B.values(), B.ld(), 
					A.values(), A.ld(), 
					tmp.values(), tmp.ld());
			C.iscale(beta);
			ops::update(T_Scalar(1), tmp, C);

		} // beta

	} else {

		throw_prop_compatibility_error(A, B, C);

	} // property combos
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Scl) \
template void mult(T_Scl, \
	op_t, const dns::XxMatrix<T_Scl>&, \
	op_t, const dns::XxMatrix<T_Scl>&, \
	T_Scl, dns::XxMatrix<T_Scl>&)
instantiate_mult(real_t);
instantiate_mult(real4_t);
instantiate_mult(complex_t);
instantiate_mult(complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Scalar>
static void trimult(T_Scalar alpha, side_t sideA, 
		op_t opA, const  dns::XxMatrix<T_Scalar>& A,
		dns::XxMatrix<T_Scalar>& B)
{
  Operation _opA(opA);

  trimat_mult_replace_check(sideA,
      A.prop(), A.nrows(), A.ncols(), _opA,
      B.prop(), B.nrows(), B.ncols());

  blas::trmm(
      static_cast<char>(sideA), A.prop().cuplo(), _opA.ctype(), 'N',
      B.nrows(), B.ncols(), alpha, A.values(), A.ld(),
      B.values(), B.ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void trimult(T_Scalar alpha, op_t opA,
		const dns::XxMatrix<T_Scalar>& A,
		dns::XxMatrix<T_Scalar>& B)
{
	trimult(alpha, side_t::Left, opA, A, B);
}
/*-------------------------------------------------*/
#define instantiate_trimult(T_Scl) \
template void trimult(T_Scl, op_t, const dns::XxMatrix<T_Scl>&, dns::XxMatrix<T_Scl>&)
instantiate_trimult(real_t);
instantiate_trimult(real4_t);
instantiate_trimult(complex_t);
instantiate_trimult(complex8_t);
#undef instantiate_trimult
/*-------------------------------------------------*/
template <typename T_Scalar>
void trimult(T_Scalar alpha,
		dns::XxMatrix<T_Scalar>& B,
		op_t opA, const dns::XxMatrix<T_Scalar>& A)
{
	trimult(alpha, side_t::Right, opA, A, B);
}
/*-------------------------------------------------*/
#define instantiate_trimult(T_Scl) \
template void trimult(T_Scl, \
		dns::XxMatrix<T_Scl>&, \
    op_t, const dns::XxMatrix<T_Scl>&)
instantiate_trimult(real_t);
instantiate_trimult(real4_t);
instantiate_trimult(complex_t);
instantiate_trimult(complex8_t);
#undef instantiate_trimult
/*-------------------------------------------------*/
template <typename T_Scalar>
static void trisol(T_Scalar alpha, side_t sideA, 
		op_t opA, const  dns::XxMatrix<T_Scalar>& A,
		dns::XxMatrix<T_Scalar>& B)
{
	Operation _opA(opA);

	trimat_mult_replace_check(sideA,
			A.prop(), A.nrows(), A.ncols(), _opA,
			B.prop(), B.nrows(), B.ncols());

	blas::trsm(
			static_cast<char>(sideA), A.prop().cuplo(), _opA.ctype(), 'N',
			B.nrows(), B.ncols(), alpha, A.values(), A.ld(),
			B.values(), B.ld());
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void trisol(T_Scalar alpha, op_t opA,
		const dns::XxMatrix<T_Scalar>& A,
		dns::XxMatrix<T_Scalar>& B)
{
	trisol(alpha, side_t::Left, opA, A, B);
}
/*-------------------------------------------------*/
#define instantiate_trisol(T_Scl) \
template void trisol(T_Scl, op_t, \
	const dns::XxMatrix<T_Scl>&, \
	dns::XxMatrix<T_Scl>&)
instantiate_trisol(real_t);
instantiate_trisol(real4_t);
instantiate_trisol(complex_t);
instantiate_trisol(complex8_t);
#undef instantiate_trisol
/*-------------------------------------------------*/
template <typename T_Scalar>
void trisol(T_Scalar alpha,
    dns::XxMatrix<T_Scalar>& B,
    op_t opA, const dns::XxMatrix<T_Scalar>& A)
{
	trisol(alpha, side_t::Right, opA, A, B);
}
/*-------------------------------------------------*/
#define instantiate_trisol(T_Scl) \
template void trisol(T_Scl, \
		dns::XxMatrix<T_Scl>&, \
    op_t, const dns::XxMatrix<T_Scl>&)
instantiate_trisol(real_t);
instantiate_trisol(real4_t);
instantiate_trisol(complex_t);
instantiate_trisol(complex8_t);
#undef instantiate_trisol
/*-------------------------------------------------*/
/*-------------------------------------------------*/
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha, op_t opA,
	const csc::XxMatrix<T_Int,T_Scalar>& A,
	const dns::XxMatrix<T_Scalar>& B,
	T_Scalar beta, dns::XxMatrix<T_Scalar>& C)
{
	if(A.prop().isSymmetric() || A.prop().isHermitian()) opA = op_t::N;

	opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);

	Operation _opA(opA);
	Operation _opB(op_t::N);

	mult_dim_check(
			A.nrows(), A.ncols(), _opA, 
			B.nrows(), B.ncols(), _opB, 
			C.nrows(), C.ncols());

	if((A.prop().isGeneral() || A.prop().isTriangular()) && B.prop().isGeneral() && C.prop().isGeneral()) {

		int_t k = (_opA.isTranspose() ? A.nrows() : A.ncols());

			blk::csc::gem_x_gem(opA, 
					C.nrows(), 
					C.ncols(), 
					k, 
					alpha,
					A.colptr(), A.rowidx(), A.values(),
					B.values(), B.ld(), 
					beta, 
					C.values(), C.ld());

	} else if(A.prop().isSymmetric() && B.prop().isGeneral() && C.prop().isGeneral()) {

		blk::csc::sym_x_gem(A.prop().uplo(),
				C.nrows(),
				C.ncols(),
				alpha,
				A.colptr(), A.rowidx(), A.values(),
				B.values(), B.ld(), 
				beta, 
				C.values(), C.ld());

	} else if(A.prop().isHermitian() && B.prop().isGeneral() && C.prop().isGeneral()) {

		blk::csc::hem_x_gem(A.prop().uplo(),
				C.nrows(),
				C.ncols(),
				alpha,
				A.colptr(), A.rowidx(), A.values(),
				B.values(), B.ld(), 
				beta, 
				C.values(), C.ld());

	} else {

		throw_prop_compatibility_error(A, B, C);

	} // property combos
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, op_t, \
	const csc::XxMatrix<T_Int,T_Scl>&, \
	const dns::XxMatrix<T_Scl>&, \
	T_Scl, dns::XxMatrix<T_Scl>&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void mult(T_Scalar alpha,
	op_t opA, const csc::XxMatrix<T_Int,T_Scalar>& A,
	op_t opB, const csc::XxMatrix<T_Int,T_Scalar>& B,
	T_Scalar beta, dns::XxMatrix<T_Scalar>& C)
{
	opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);
	opB = (TypeTraits<T_Scalar>::is_real() && opB == op_t::C ? op_t::T : opB);

	Operation _opA(opA);
	Operation _opB(opB);

	mult_dim_check(
			A.nrows(), A.ncols(), _opA, 
			B.nrows(), B.ncols(), _opB, 
			C.nrows(), C.ncols());

	if(A.prop().isGeneral() && A.prop().isGeneral() && A.prop().isGeneral()) {

		int_t k = (_opA.isTranspose() ? A.nrows() : A.ncols());

		blk::csc::gem_x_gem(C.nrows(), C.ncols(), k, alpha,
				opA, A.colptr(), A.rowidx(), A.values(),
				opB, B.colptr(), B.rowidx(), B.values(),
				beta, C.values(), C.ld());

	} else {

		throw_prop_compatibility_error(A, B, C);

	} // property combos
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template void mult(T_Scl, \
		op_t, const csc::XxMatrix<T_Int,T_Scl>&, \
		op_t, const csc::XxMatrix<T_Int,T_Scl>&, \
		T_Scl, dns::XxMatrix<T_Scl>&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
csc::XxMatrix<T_Int,T_Scalar> mult(T_Scalar alpha,
    op_t opA, const csc::XxMatrix<T_Int,T_Scalar>& A,
    op_t opB, const csc::XxMatrix<T_Int,T_Scalar>& B)
{
	opA = (TypeTraits<T_Scalar>::is_real() && opA == op_t::C ? op_t::T : opA);
	opB = (TypeTraits<T_Scalar>::is_real() && opB == op_t::C ? op_t::T : opB);

	Operation _opA(opA);
	Operation _opB(opB);

	int_t m = (_opA.isTranspose() ? A.ncols() : A.nrows());
	int_t n = (_opB.isTranspose() ? B.nrows() : B.ncols());
	int_t k = (_opA.isTranspose() ? A.nrows() : A.ncols());

	mult_dim_check(
			A.nrows(), A.ncols(), _opA, 
			B.nrows(), B.ncols(), _opB, 
			m, n);

	csc::XxMatrix<T_Int,T_Scalar> ret;

	if(A.prop().isGeneral() && B.prop().isGeneral()) {

		int_t    *colptrC = nullptr;
		int_t    *rowidxC = nullptr;
		T_Scalar *valuesC = nullptr;

		blk::csc::gem_x_gem(m, n, k,
				opA, A.colptr(), A.rowidx(), A.values(),
				opB, B.colptr(), B.rowidx(), B.values(),
				&colptrC, &rowidxC, &valuesC);

		ret = csc::XxMatrix<T_Int,T_Scalar>(m, n, colptrC, rowidxC, valuesC, true);
		ret.iscale(alpha);

	} else {

		throw_prop_compatibility_error(A, B);

	} // property combos

	return ret;
}
/*-------------------------------------------------*/
#define instantiate_mult(T_Int, T_Scl) \
template csc::XxMatrix<T_Int,T_Scl> mult(T_Scl, \
    op_t, const csc::XxMatrix<T_Int,T_Scl>&, \
    op_t, const csc::XxMatrix<T_Int,T_Scl>&)
instantiate_mult(int_t, real_t);
instantiate_mult(int_t, real4_t);
instantiate_mult(int_t, complex_t);
instantiate_mult(int_t, complex8_t);
#undef instantiate_mult
/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/
