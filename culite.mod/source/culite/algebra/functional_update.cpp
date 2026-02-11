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
#include "culite/algebra/functional_update.hpp"

// system

// 3rd
#include <cla3p/checks/basic_checks.hpp>

// culite
#include "culite/error/exceptions.hpp"
#include "culite/dense/dns_xxvector.hpp"
#include "culite/dense/dns_xxmatrix.hpp"
#include "culite/sparse/csr_xxmatrix.hpp"
#include "culite/sparse/csc_xxmatrix.hpp"

/*-------------------------------------------------*/
namespace culite {
namespace ops {
/*-------------------------------------------------*/
template <typename T_Scalar>
void update(T_Scalar alpha, const dns::XxVector<T_Scalar>& x, dns::XxVector<T_Scalar>& y, CuBlasHandler& cublasHandler)
{
	::cla3p::similarity_dim_check(x.size(), y.size());

	cublasHandler.axpy(x.size(), &alpha, x.values(), 1, y.values(), 1);
}
/*-------------------------------------------------*/
#define instantiate_update(T_Scl) \
template void update(T_Scl, const dns::XxVector<T_Scl>&, dns::XxVector<T_Scl>&, CuBlasHandler&)
instantiate_update(real_t);
instantiate_update(real4_t);
instantiate_update(complex_t);
instantiate_update(complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
template <typename T_Scalar>
void update(::cla3p::op_t opA,
            T_Scalar alpha,
            const dns::XxMatrix<T_Scalar>& A,
            dns::XxMatrix<T_Scalar>& B,
            CuBlasHandler& cublasHandler)
{
	::cla3p::similarity_check(opA, A, ::cla3p::op_t::N, B);

    T_Scalar beta = makeScalar<T_Scalar>(1);

    if(B.prop().isGeneral()) {

        cublasHandler.geam<T_Scalar>(opA,
				                     ::cla3p::op_t::N,
							         B.nrows(), B.ncols(),
				                     &alpha, A.values(), A.ld(),
				                     &beta, B.values(), B.ld(),
				                     B.values(), B.ld());

    } else {

        throw err::CudaException("Unsupported matrix property for updating. Only General (dense) matrices are currently supported.");

    } // prop
}
/*-------------------------------------------------*/
#define instantiate_update(T_Scl) \
template void update(::cla3p::op_t, T_Scl, const dns::XxMatrix<T_Scl>&, dns::XxMatrix<T_Scl>&, CuBlasHandler&)
instantiate_update(real_t);
instantiate_update(real4_t);
instantiate_update(complex_t);
instantiate_update(complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void update(T_Scalar /*alpha*/,
            const csr::XxMatrix<T_Int,T_Scalar>& /*A*/,
            csr::XxMatrix<T_Int,T_Scalar>& /*B*/,
            CuSparseHandler& /*cuSparseHandler*/)
{
    throw err::CudaException("CSR matrix update not supported yet.");
}
/*-------------------------------------------------*/
#define instantiate_update(T_Int,T_Scl) \
template void update(T_Scl, const csr::XxMatrix<T_Int,T_Scl>&, csr::XxMatrix<T_Int,T_Scl>&, CuSparseHandler&)
instantiate_update(int_t,real_t);
instantiate_update(int_t,real4_t);
instantiate_update(int_t,complex_t);
instantiate_update(int_t,complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
template <typename T_Int, typename T_Scalar>
void update(T_Scalar /*alpha*/,
            const csc::XxMatrix<T_Int,T_Scalar>& /*A*/,
            csc::XxMatrix<T_Int,T_Scalar>& /*B*/,
            CuSparseHandler& /*cuSparseHandler*/)
{
    throw err::CudaException("CSC matrix update not supported yet.");
}
/*-------------------------------------------------*/
#define instantiate_update(T_Int,T_Scl) \
template void update(T_Scl, const csc::XxMatrix<T_Int,T_Scl>&, csc::XxMatrix<T_Int,T_Scl>&, CuSparseHandler&)
instantiate_update(int_t,real_t);
instantiate_update(int_t,real4_t);
instantiate_update(int_t,complex_t);
instantiate_update(int_t,complex8_t);
#undef instantiate_update
/*-------------------------------------------------*/
} // namespace ops
} // namespace culite
/*-------------------------------------------------*/
