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
#include "cla3p/algebra/functional_outer.hpp"

// system

// 3rd

// cla3p
#include "cla3p/error/exceptions.hpp"
#include "cla3p/checks/outer_checks.hpp"
#include "cla3p/checks/hermitian_coeff_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
static void outerx(bool conjop, T_Scalar alpha, 
		const dns::XxVector<T_Scalar>& X, 
		const dns::XxVector<T_Scalar>& Y, 
		dns::XxMatrix<T_Scalar>& A)
{
	conjop = (TypeTraits<T_Scalar>::is_real() ? false : conjop);

	outer_product_consistency_check(conjop, A.nrows(), A.ncols(), A.prop(), X.size(), Y.size());
	hermitian_coeff_check(A.prop(), alpha);

	if(A.prop().isGeneral()) {

		if(conjop) blas::gerc(X.size(), Y.size(), alpha, X.values(), 1, Y.values(), 1, A.values(), A.ld());
		else       blas::ger (X.size(), Y.size(), alpha, X.values(), 1, Y.values(), 1, A.values(), A.ld());

	} else if(A.prop().isSymmetric()) {

		if(X.values() == Y.values()) {
			blas::syr(A.prop().cuplo(), X.size(), alpha, X.values(), 1, A.values(), A.ld());
		} else {
			blas::gemmt(A.prop().cuplo(), 'N', 'T', X.size(), 1, alpha, X.values(), X.size(), Y.values(), Y.size(), 1, A.values(), A.ld());
		}

	} else if(A.prop().isHermitian()) {

		if(X.values() == Y.values()) {
			blas::her(A.prop().cuplo(), X.size(), arith::getRe(alpha), X.values(), 1, A.values(), A.ld());
		} else {
			blas::gemmt(A.prop().cuplo(), 'N', 'C', X.size(), 1, alpha, X.values(), X.size(), Y.values(), Y.size(), 1, A.values(), A.ld());
		}

	} else {

		throw err::Exception();

	} // valid props
}
/*-------------------------------------------------*/
template <typename T_Scalar>
void outer(T_Scalar alpha, 
		const dns::XxVector<T_Scalar>& X, 
		const dns::XxVector<T_Scalar>& Y, 
		dns::XxMatrix<T_Scalar>& A)
{
	outerx(false, alpha, X, Y, A);
}
/*-------------------------------------------------*/
#define instantiate_outer(T_Scl) \
template void outer(T_Scl, const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&, dns::XxMatrix<T_Scl>&)
instantiate_outer(real_t);
instantiate_outer(real4_t);
instantiate_outer(complex_t);
instantiate_outer(complex8_t);
#undef instantiate_outer
/*-------------------------------------------------*/
template <typename T_Scalar>
void outerc(T_Scalar alpha, 
		const dns::XxVector<T_Scalar>& X, 
		const dns::XxVector<T_Scalar>& Y, 
		dns::XxMatrix<T_Scalar>& A)
{
	outerx(true, alpha, X, Y, A);
}
/*-------------------------------------------------*/
#define instantiate_outerc(T_Scl) \
template void outerc(T_Scl, const dns::XxVector<T_Scl>&, const dns::XxVector<T_Scl>&, dns::XxMatrix<T_Scl>&)
instantiate_outerc(real_t);
instantiate_outerc(real4_t);
instantiate_outerc(complex_t);
instantiate_outerc(complex8_t);
#undef instantiate_outerc
/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/
