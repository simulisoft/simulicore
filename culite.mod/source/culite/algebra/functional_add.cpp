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
#include "culite/algebra/functional_add.hpp"

// system

// 3rd
#include <cla3p/types/property.hpp>
#include <cla3p/checks/basic_checks.hpp>

// culite
#include "culite/error/exceptions.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
dns::XxVector<T_Scalar> add(
		T_Scalar alpha, const dns::XxVector<T_Scalar>& x,
		T_Scalar beta , const dns::XxVector<T_Scalar>& y,
        CuBlasHandler& cublasHandler)
{
	::cla3p::similarity_dim_check(x.size(), y.size());

	dns::XxVector<T_Scalar> ret(x.size());

    cublasHandler.geam(::cla3p::op_t::N,
                       ::cla3p::op_t::N,
                       ret.size(), 1,
                       &alpha, x.values(), x.size(),
                       &beta,  y.values(), y.size(),
                       ret.values(), ret.size());

	return ret;
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template dns::XxVector<T_Scl> add( \
		T_Scl, const dns::XxVector<T_Scl>&, \
		T_Scl, const dns::XxVector<T_Scl>&, \
        CuBlasHandler&)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
template <typename T_Scalar>
 dns::XxMatrix<T_Scalar> add(
		::cla3p::op_t opA, T_Scalar alpha, const dns::XxMatrix<T_Scalar>& A,
		::cla3p::op_t opB, T_Scalar beta , const dns::XxMatrix<T_Scalar>& B,
        CuBlasHandler& cublasHandler)
{
	::cla3p::similarity_check(opA, A, opB, B);

	int_t m = A.nrows();
	int_t n = A.ncols();
	::cla3p::Property pr = A.prop();

    dns::XxMatrix<T_Scalar> ret(m, n, pr);

    if(pr.isGeneral()) {

        cublasHandler.geam(opA,
                           opB,
                           m, n,
                           &alpha, A.values(), A.ld(),
                           &beta,  B.values(), B.ld(),
                           ret.values(), ret.ld());

    } else {

        throw err::CudaException("Unsupported matrix property for addition. Only General (dense) matrices are currently supported.");

    } // prop

	return ret;
}
/*-------------------------------------------------*/
#define instantiate_add(T_Scl) \
template dns::XxMatrix<T_Scl> add( \
		::cla3p::op_t, T_Scl, const dns::XxMatrix<T_Scl>&, \
		::cla3p::op_t, T_Scl, const dns::XxMatrix<T_Scl>&, \
        CuBlasHandler&)
instantiate_add(real_t);
instantiate_add(real4_t);
instantiate_add(complex_t);
instantiate_add(complex8_t);
#undef instantiate_add
/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/
