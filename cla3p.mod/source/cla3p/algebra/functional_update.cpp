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
#include "cla3p/algebra/functional_update.hpp"

// system

// 3rd

// cla3p
#include "cla3p/checks/basic_checks.hpp"
#include "cla3p/proxies/blas_proxy.hpp"
#include "cla3p/bulk/dns_math.hpp"
#include "cla3p/dense/dns_xxvector.hpp"
#include "cla3p/dense/dns_xxmatrix.hpp"
#include "cla3p/sparse/csc_xxmatrix.hpp"
#include "cla3p/algebra/functional_add.hpp"

/*-------------------------------------------------*/
namespace cla3p {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
void update(T_Scalar alpha, const dns::XxVector<T_Scalar>& X, dns::XxVector<T_Scalar>& Y)
{
	similarity_dim_check(X.size(), Y.size());
	blas::axpy(X.size(), alpha, X.values(), 1, Y.values(), 1);
}
/*-------------------------------------------------*/
#define instantiate_update(T_Scl) \
template void update(T_Scl, const dns::XxVector<T_Scl>&, dns::XxVector<T_Scl>&)
instantiate_update(real_t);
instantiate_update(real4_t);
instantiate_update(complex_t);
instantiate_update(complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
template <typename T_Scalar>
void update(T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A, dns::XxMatrix<T_Scalar>& B)
{
	similarity_check(
			A.prop(), A.nrows(), A.ncols(),
			B.prop(), B.nrows(), B.ncols());

	blk::dns::update(A.prop().uplo(), A.nrows(), A.ncols(), alpha, 
			A.values(), A.ld(),
			B.values(), B.ld());
}
/*-------------------------------------------------*/
#define instantiate_update(T_Scl) \
template void update(T_Scl, const dns::XxMatrix<T_Scl>&, dns::XxMatrix<T_Scl>&)
instantiate_update(real_t);
instantiate_update(real4_t);
instantiate_update(complex_t);
instantiate_update(complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void update(T_Scalar alpha, const csc::XxMatrix<T_Int,T_Scalar>& A, csc::XxMatrix<T_Int,T_Scalar>& B)
{
	similarity_check(
			A.prop(), A.nrows(), A.ncols(),
			B.prop(), B.nrows(), B.ncols());

	csc::XxMatrix<T_Int,T_Scalar> tmp = add(alpha, A, T_Scalar(1), B);
	B.clear();
	B = tmp.move();
}
/*-------------------------------------------------*/
#define instantiate_update(T_Int,T_Scl) \
template void update(T_Scl, const csc::XxMatrix<T_Int,T_Scl>&, csc::XxMatrix<T_Int,T_Scl>&)
instantiate_update(int_t, real_t);
instantiate_update(int_t, real4_t);
instantiate_update(int_t, complex_t);
instantiate_update(int_t, complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
} // namespace ops
} // namespace cla3p
/*-------------------------------------------------*/
